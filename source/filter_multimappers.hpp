#ifndef FILTER_MULTIMAPPERS_H
#define FILTER_MULTIMAPPERS_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_multimappers(chimeric_alignments_t& chimeric_alignments, fusions_t& fusions, const exon_annotation_index_t& exon_annotation_index, const assembly_t& assembly,disjoint_set& homolog_union);

#endif /* FILTER_MULTIMAPPERS_H */
