// Read-only stock-release observation. Each wrapper forwards to its real
// function exactly once and records the return or emitted value afterward.
#include <JetJvtEfficiency/NNJvtSelectionTool.h>
#include <METUtilities/METMaker.h>
#include <METUtilities/METSignificance.h>
#include <METUtilities/METSystematicsTool.h>
#include <MetAnalysisAlgorithms/MetBuilderAlg.h>
#include <MetAnalysisAlgorithms/MetMakerAlg.h>
#include <MetAnalysisAlgorithms/MetSignificanceAlg.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <xAODMissingET/MissingETAssociationHelper.h>

#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>

namespace {
thread_local long event_entry = -1;
thread_local std::size_t maker_state = 0;
thread_local std::size_t significance_state = 0;
thread_local bool in_maker = false;
thread_local bool in_significance = false;
thread_local bool in_correction = false;
thread_local bool in_seed = false;

template <class Signature>
Signature release_function(const char* symbol) {
  auto* address = dlsym(RTLD_NEXT, symbol);
  if (!address) {
    const char* stock_libraries[] = {
        "libMetAnalysisAlgorithmsLib.so", "libMETUtilitiesLib.so",
        "libPATInterfaces.so", "libMathCore.so",
        "libJetJvtEfficiencyLib.so", "libxAODMissingET.so"};
    for (const char* library : stock_libraries) {
      auto* handle = dlopen(library, RTLD_LAZY | RTLD_NOLOAD);
      if (!handle) {
        continue;
      }
      address = dlsym(handle, symbol);
      dlclose(handle);
      if (address) {
        break;
      }
    }
  }
  if (!address) {
    std::fprintf(stderr, "MET stock trace could not resolve %s: %s\n", symbol, dlerror());
    std::abort();
  }
  return reinterpret_cast<Signature>(address);
}

std::FILE* trace_file() {
  static std::FILE* output = []() {
    const char* path = std::getenv("MET_STOCK_TRACE");
    auto* file = path ? std::fopen(path, "wx") : nullptr;
    if (!file) {
      std::fputs("MET stock trace requires a fresh MET_STOCK_TRACE file\n", stderr);
      std::abort();
    }
    return file;
  }();
  return output;
}

template <class... Args>
void trace(std::size_t state, const char* action, const char* format, Args... fields) {
  auto* output = trace_file();
  std::fprintf(output, "%ld\t%zu\t%s\t", event_entry, state, action);
  std::fprintf(output, format, fields...);
  std::fputc('\n', output);
}
}  // namespace

StatusCode CP::MetMakerAlg::execute(const EventContext& context) {
  using Function = StatusCode (*)(MetMakerAlg*, const EventContext&);
  static auto original = release_function<Function>("_ZN2CP11MetMakerAlg7executeERK12EventContext");
  ++event_entry;
  maker_state = 0;
  in_maker = true;
  const auto status = original(this, context);
  in_maker = false;
  trace(maker_state, "makerAlg", "%lu", static_cast<unsigned long>(status));
  return status;
}

StatusCode CP::MetBuilderAlg::execute(const EventContext& context) {
  using Function = StatusCode (*)(MetBuilderAlg*, const EventContext&);
  static auto original = release_function<Function>("_ZN2CP13MetBuilderAlg7executeERK12EventContext");
  const auto status = original(this, context);
  trace(maker_state, "builderAlg", "%lu", static_cast<unsigned long>(status));
  return status;
}

StatusCode CP::MetSignificanceAlg::execute(const EventContext& context) {
  using Function = StatusCode (*)(MetSignificanceAlg*, const EventContext&);
  static auto original = release_function<Function>("_ZN2CP18MetSignificanceAlg7executeERK12EventContext");
  significance_state = 0;
  in_significance = true;
  const auto status = original(this, context);
  in_significance = false;
  trace(significance_state, "significanceAlg", "%lu", static_cast<unsigned long>(status));
  return status;
}

StatusCode met::METMaker::rebuildMET(
    const std::string& name, xAOD::Type::ObjectType type, xAOD::MissingETContainer* met,
    const xAOD::IParticleContainer* selected, xAOD::MissingETAssociationHelper& helper,
    MissingETBase::UsageHandler::Policy scale) const {
  using Function = StatusCode (*)(
      const METMaker*, const std::string&, xAOD::Type::ObjectType, xAOD::MissingETContainer*,
      const xAOD::IParticleContainer*, xAOD::MissingETAssociationHelper&, MissingETBase::UsageHandler::Policy);
  static auto original = release_function<Function>(
      "_ZNK3met8METMaker10rebuildMETERKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEN8xAODType10ObjectTypeEPN4xAOD21MissingETContainer_v1EPK10DataVectorINSB_9IParticleEN16DataModel_detail6NoBaseEERNSB_26MissingETAssociationHelperEN13MissingETBase12UsageHandler6PolicyE");
  const auto status = original(this, name, type, met, selected, helper, scale);
  if (in_maker) {
    trace(maker_state, "rebuildMET", "%s\t%zu\t%lu",
          name.c_str(), selected->size(), static_cast<unsigned long>(status));
  }
  return status;
}

StatusCode met::METMaker::rebuildJetMET(
    const std::string& jet, const std::string& soft, xAOD::MissingETContainer* met,
    const xAOD::JetContainer* jets, const xAOD::MissingETContainer* core,
    xAOD::MissingETAssociationHelper& helper, bool do_jvt) const {
  using Function = StatusCode (*)(
      const METMaker*, const std::string&, const std::string&, xAOD::MissingETContainer*,
      const xAOD::JetContainer*, const xAOD::MissingETContainer*, xAOD::MissingETAssociationHelper&, bool);
  static auto original = release_function<Function>(
      "_ZNK3met8METMaker13rebuildJetMETERKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_PN4xAOD21MissingETContainer_v1EPK10DataVectorINS9_6Jet_v1ESC_INS9_9IParticleEN16DataModel_detail6NoBaseEEEPKSA_RNS9_26MissingETAssociationHelperEb");
  const auto status = original(this, jet, soft, met, jets, core, helper, do_jvt);
  if (in_maker) {
    trace(maker_state, "rebuildJetMET", "%s\t%s\t%zu\t%d\t%lu",
          jet.c_str(), soft.c_str(), jets->size(), do_jvt, static_cast<unsigned long>(status));
  }
  return status;
}

StatusCode CP::SystematicsTool::applySystematicVariation(const CP::SystematicSet& state) {
  using Function = StatusCode (*)(SystematicsTool*, const CP::SystematicSet&);
  static auto original = release_function<Function>("_ZN2CP15SystematicsTool24applySystematicVariationERKNS_13SystematicSetE");
  const auto status = original(this, state);
  if (in_maker) {
    trace(maker_state, "applySystematicVariation", "%s\t%lu", state.name().c_str(), static_cast<unsigned long>(status));
  }
  return status;
}

void met::METSystematicsTool::setRandomSeed(unsigned long seed) const {
  using Function = void (*)(const METSystematicsTool*, unsigned long);
  static auto original = release_function<Function>("_ZNK3met18METSystematicsTool13setRandomSeedEm");
  in_seed = in_maker;
  original(this, seed);
  in_seed = false;
  if (in_maker) {
    trace(maker_state, "setRandomSeed", "%lu", seed);
  }
}

CP::CorrectionCode met::METSystematicsTool::applyCorrection(
    xAOD::MissingET& soft, const xAOD::MissingETAssociationHelper& helper) const {
  using Function = CP::CorrectionCode (*)(
      const METSystematicsTool*, xAOD::MissingET&, const xAOD::MissingETAssociationHelper&);
  static auto original = release_function<Function>(
      "_ZNK3met18METSystematicsTool15applyCorrectionERN4xAOD12MissingET_v1ERKNS1_26MissingETAssociationHelperE");
  if (!in_maker) {
    return original(this, soft, helper);
  }
  const double input_x = soft.mpx();
  const double input_y = soft.mpy();
  const double input_sumet = soft.sumet();
  in_correction = true;
  const auto status = original(this, soft, helper);
  in_correction = false;
  trace(maker_state, "applyCorrection", "%d\t%a\t%a\t%a\t%a\t%a\t%a",
        static_cast<int>(status.code()), input_x, input_y, input_sumet,
        static_cast<double>(soft.mpx()), static_cast<double>(soft.mpy()), static_cast<double>(soft.sumet()));
  ++maker_state;
  return status;
}

StatusCode met::METSignificance::varianceMET(
    xAOD::MissingETContainer* met, float average_mu, const std::string& jet,
    const std::string& soft, const std::string& total) {
  using Function = StatusCode (*)(METSignificance*, xAOD::MissingETContainer*, float,
                                 const std::string&, const std::string&, const std::string&);
  static auto original = release_function<Function>(
      "_ZN3met15METSignificance11varianceMETEPN4xAOD21MissingETContainer_v1EfRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESB_SB_");
  const auto status = original(this, met, average_mu, jet, soft, total);
  if (in_significance) {
    trace(significance_state, "varianceMET", "%lu\t%a",
          static_cast<unsigned long>(status), static_cast<double>(average_mu));
    ++significance_state;
  }
  return status;
}

double TRandom::Gaus(double mean, double standard_deviation) {
  using Function = double (*)(TRandom*, double, double);
  static auto original = release_function<Function>("_ZN7TRandom4GausEdd");
  const double draw = original(this, mean, standard_deviation);
  if (in_correction) {
    trace(maker_state, "gaussian", "%a\t%a\t%a", mean, standard_deviation, draw);
  }
  return draw;
}

bool CP::NNJvtSelectionTool::select(const xAOD::IParticle* jet) const {
  using Function = bool (*)(const NNJvtSelectionTool*, const xAOD::IParticle*);
  static auto original = release_function<Function>("_ZNK2CP18NNJvtSelectionTool6selectEPKN4xAOD9IParticleE");
  const bool accept = original(this, jet);
  if (in_maker) {
    trace(maker_state, "jetJvt", "%zu\t%d", jet->index(), accept);
  }
  return accept;
}

void xAOD::MissingETAssociationHelper::setObjSelectionFlag(
    const MissingETAssociation_v1* association, std::size_t index, bool selected) {
  using Function = void (*)(MissingETAssociationHelper*, const MissingETAssociation_v1*, std::size_t, bool);
  static auto original = release_function<Function>("_ZN4xAOD26MissingETAssociationHelper19setObjSelectionFlagEPKNS_23MissingETAssociation_v1Emb");
  original(this, association, index, selected);
  if (in_maker) {
    trace(maker_state, "associationFlag", "%zu\t%zu\t%d", association->index(), index, selected);
  }
}

void TRandom3::SetSeed(unsigned long seed) {
  using Function = void (*)(TRandom3*, unsigned long);
  static auto original = release_function<Function>("_ZN8TRandom37SetSeedEm");
  original(this, seed);
  if (in_seed) {
    trace(maker_state, "rootSeed", "%lu", seed);
  }
}
