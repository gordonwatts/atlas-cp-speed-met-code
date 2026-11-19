// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_COLUMNAR_TRANSPORT_H
#define ATLAS_MET_COLUMNAR_TRANSPORT_H

#include <stddef.h>
#include <stdint.h>

extern "C" {
struct MetColumn {
    const char* name;
    const char* dtype;
    size_t size;
    const void* data;
};
struct MetColumnar;
const char* met_columnar_error();
MetColumnar* met_columnar_open(const char* provider, const char* maker, const char* nnjvt);
void met_columnar_close(MetColumnar*);
int met_columnar_configure(MetColumnar*, const char* soft, const char* significance,
                          const char* preparation_provider, const char* preparation);
void met_columnar_reset(MetColumnar*);
int met_prep_evaluate(const MetColumnar*, const char* operation, size_t rows,
                      const MetColumn*, size_t columns, double* output, size_t components);
int met_columnar_full(MetColumnar*, size_t events, uint64_t first_entry,
                      const char* const* states, size_t state_count,
                      const MetColumn*, size_t columns, float* terms, uint64_t* sources,
                      float* significance, const char* trace_file);

int met_decode_jet_links(const uint8_t* raw, size_t raw_size,
                         const int64_t* offsets, size_t events, int64_t* decoded_offsets,
                         uint32_t* key, uint32_t* index, size_t capacity);
int met_columnar_raw_smoke(const MetColumnar*, size_t events, const MetColumn*, size_t columns,
                           float* terms, size_t terms_count, uint64_t* association_flags,
                           size_t association_count, float* jet_scores, int8_t* jet_jvt,
                           size_t jet_count);
}
#endif
