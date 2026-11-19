// Small C ABI for a release-free Python diagnostic of the shared evaluator.
#include <correction.h>
#include <correction_plugin.h>
#include "small_jvt.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string error_message;

struct Reader {
  correction::PluginRegistry registry;
  std::unique_ptr<correction::CorrectionSet> set;

  explicit Reader(const char* payload) {
    register_small_jvt(registry);
    set = correction::CorrectionSet::from_file(payload, registry);
  }
};
}  // namespace

extern "C" const char* met_jvt_error() {
  return error_message.c_str();
}

extern "C" Reader* met_jvt_open(const char* payload) {
  try {
    return new Reader(payload);
  } catch (const std::exception& error) {
    error_message = error.what();
    return nullptr;
  }
}

extern "C" void met_jvt_close(Reader* handle) {
  delete handle;
}

extern "C" int met_jvt_evaluate(
    const Reader* handle, std::size_t size, const double* pt, const double* eta,
    const double* sum_pt_trk, const double* jvf_corr, const std::int64_t* vertex,
    const std::int64_t* feature_state, double* output) {
  try {
    const auto update = handle->set->at("small_nnjvt_update");
    const auto select = handle->set->at("small_nnjvt_select");
    std::array<double, 4> update_result{};
    std::array<double, 2> select_result{};
    std::vector<correction::Variable::Type> update_inputs(6), select_inputs(4);
    for (std::size_t row = 0; row < size; ++row) {
      update_inputs[0] = pt[row];
      update_inputs[1] = eta[row];
      update_inputs[2] = sum_pt_trk[row];
      update_inputs[3] = jvf_corr[row];
      update_inputs[4] = vertex[row];
      update_inputs[5] = feature_state[row];
      update->evaluate_array(update_inputs, update_result.data(), update_result.size());
      if (update_result[3] != 1.0) {
        throw std::runtime_error("Release diagnostic saw skipped/unavailable NNJVT input on row " + std::to_string(row));
      }
      select_inputs[0] = pt[row];
      select_inputs[1] = eta[row];
      select_inputs[2] = update_result[0];
      select_inputs[3] = std::int64_t{1};
      select->evaluate_array(select_inputs, select_result.data(), select_result.size());
      if (select_result[1] != 1.0) {
        throw std::runtime_error("Release diagnostic saw unexpected selector status");
      }
      output[3 * row] = update_result[0];
      output[3 * row + 1] = update_result[2];
      output[3 * row + 2] = select_result[0];
    }
    return 0;
  } catch (const std::exception& error) {
    error_message = error.what();
    return 1;
  }
}
