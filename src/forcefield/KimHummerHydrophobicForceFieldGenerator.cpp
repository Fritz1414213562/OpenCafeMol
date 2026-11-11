#include "KimHummerHydrophobicForceFieldGenerator.hpp"

#include "ForceFieldIDGenerator.hpp"
#include "InteractionGroup.hpp"
#include "src/util/Logger.hpp"
#include "src/util/Constants.hpp"

#include <iostream>
#include <unordered_map>

KimHummerHydrophobicForceFieldGenerator::KimHummerHydrophobicForceFieldGenerator(
    const double cutoff_ratio,
    const double epsilon,
	const double alpha,
	const double temperature,
    const std::vector<std::optional<double>> sigmas,
	const std::vector<std::optional<std::string>> restypes,
    const index_pairs_type& ignore_list, const bool use_periodic,
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
    const std::vector<std::optional<std::string>> group_vec)
    : cutoff_ratio_(cutoff_ratio), epsilon_(epsilon), alpha_(alpha),
	  RT_(temperature * Constant::Na * Constant::kB * 1e-3),
      sigmas_(sigmas), restypes_(restypes), ignore_list_(ignore_list),
      use_periodic_(use_periodic), ffgen_id_(fmt::format("KMHPS{}", ffid.gen()))
{

    log_assert(this->sigmas_.size() == this->restypes_.size(),
        "LennardJonesAttractive: sigma.size(={}) must be the same as restype.size(={})",
        sigmas_.size(), restypes_.size());

    for(std::size_t i=0; i<sigmas_.size(); ++i)
    {
        log_assert(sigmas_.at(i).has_value() == restypes_.at(i).has_value(),
            "KimHummerHydrophobic: if particle i(={}) has sigma, "
            "it should have restype too, and vice versa.", i);
    }

    interaction_groups_ = extract_interaction_group(
            this->restypes_, ignore_group_pairs, group_vec);

    cutoff_correction_ =
        4.0*(std::pow(1.0 / cutoff_ratio, 12) - std::pow(1.0 / cutoff_ratio, 6));
}

std::unique_ptr<OpenMM::Force> KimHummerHydrophobicForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
		"step(threshold - r) * (LJ + (1 - lambda) * epsilon) + step(r - threshold) * lambda * LJ;"
		"LJ         = 4 * epsilon * (sigma_r_12 - sigma_r_6);"
        "sigma_r_12 = sigma_r_6^2;"
        "sigma_r_6  = sigma_r^6;"
        "sigma_r    = sigma/r;"
        "threshold  = sigma*2^(1/6);"
        "sigma      = ({id}_sigma1 + {id}_sigma2) * 0.5;"
		"lambda     = 2 * step({id}_epsilon - mj_eps) - 1;"
        "epsilon    = abs({id}_alpha * (mj_eps - {id}_epsilon));"
		"mj_eps     = {id}_table({id}_restype1, {id}_restype2);"
		,
        fmt::arg("id", ffgen_id_));
    auto kh_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);
	auto mj_matrix = std::make_unique<OpenMM::Discrete2DFunction>(n_restype, n_restype, MiyazawaJerniganContactPotentialParameter::epsilon_MJ(this->RT_));
	auto mj_type2idx = MiyazawaJerniganContactPotentialParameter::type2idx();

	kh_ff->addTabulatedFunction(fmt::format("{}_table",      ffgen_id_), mj_matrix.release());
    kh_ff->addGlobalParameter(fmt::format("{}_epsilon",      ffgen_id_), epsilon_);
    kh_ff->addGlobalParameter(fmt::format("{}_alpha",        ffgen_id_), alpha_);
    kh_ff->addPerParticleParameter(fmt::format("{}_sigma",   ffgen_id_));
	kh_ff->addPerParticleParameter(fmt::format("{}_restype", ffgen_id_));

    double max_sigma        = std::numeric_limits<double>::min();
    double second_max_sigma = std::numeric_limits<double>::min();
    for(std::size_t idx=0; idx<sigmas_.size(); ++idx)
    {
        const std::optional<double>& sigma   = sigmas_[idx];
		const std::optional<std::string>& restype = restypes_[idx];
        if (sigma && restype)
        {
            double sigma_val = sigma.value();
			std::size_t res_idx = mj_type2idx[restype.value()];
            kh_ff->addParticle({sigma_val, res_idx});

            if(max_sigma <= sigma_val)
            {
                second_max_sigma = max_sigma;
                max_sigma        = sigma_val;
            }
            else if(second_max_sigma <= sigma_val)
            {
                second_max_sigma = sigma_val;
            }
        }
        else if(!sigma && !restype)
        {
            kh_ff->addParticle({std::numeric_limits<double>::quiet_NaN(),
								std::numeric_limits<double>::quiet_NaN()});
        }
        else
        {
            throw std::runtime_error(
                "[error] KimHummerHydrophobicForceFieldGenerator : "
                "The parameter set is epsilon, sigma and restype. incomplete parameter set "
                "was given for particle idx " + std::to_string(idx) + ".");
        }
    }

    // if interaction_group size is 0, no interaction group will be added,
    // so all the particle in the system will be considered as participant
    for(const auto& group_pair : interaction_groups_)
    {
        kh_ff->addInteractionGroup(group_pair.first, group_pair.second);
    }

    // set pbc condition
    if(use_periodic_)
    {
        kh_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
    }
    else
    {
        kh_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }

    // set cutoff
    const double cutoff_distance =
        (max_sigma + second_max_sigma) * 0.5 * cutoff_ratio_;
    std::cerr << "    KimHummerHydrophobic        : cutoff distance is "
              << cutoff_distance << " nm" << std::endl;
    kh_ff->setCutoffDistance(cutoff_distance);
    kh_ff->setUseSwitchingFunction(true);
    kh_ff->setSwitchingDistance(0.9*cutoff_distance);

    // set exclusion list
    for(const auto& pair : ignore_list_)
    {
        kh_ff->addExclusion(pair.first, pair.second);
    }

    return kh_ff;
}
