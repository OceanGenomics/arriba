#include <climits>
#include <unordered_map>
#include <unordered_set>
#include "common.hpp"
#include "annotation.hpp"
#include "assembly.hpp"
#include "filter_multimappers.hpp"
#include <iostream>
#include <random>

using namespace std;


bool is_gap_at_splice_site(const position_t position, const direction_t direction, const gene_set_t& genes, const exon_annotation_index_t& exon_annotation_index) {
	for (gene_set_t::const_iterator gene = genes.begin(); gene != genes.end(); ++gene)
		if (is_breakpoint_spliced(*gene, direction, position, exon_annotation_index))
			return true;
	return false;
}

int calculate_segment_score(const alignment_t& alignment, const string& sequence, const exon_annotation_index_t& exon_annotation_index, const assembly_t& assembly) {

	if (assembly.find(alignment.contig) == assembly.end())
		return 0;

	int score = 0;
	position_t reference_position = alignment.start;
	position_t read_position = 0;
	for (unsigned int i = 0; i < alignment.cigar.size(); ++i) {
		switch (alignment.cigar.operation(i)) {
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP: // there is no difference between soft- and hard-clipping, because we take the sequence from the SPLIT_READ, which is never hard-clipped
				read_position += alignment.cigar.op_length(i);
				break;
			case BAM_CDEL:
				score--;
				reference_position += alignment.cigar.op_length(i);
				break;
			case BAM_CREF_SKIP:
				if (!is_gap_at_splice_site(reference_position, DOWNSTREAM, alignment.genes, exon_annotation_index) ||
				    !is_gap_at_splice_site(reference_position + alignment.cigar.op_length(i), UPSTREAM, alignment.genes, exon_annotation_index))
					score--; // penalize reference skips except at splice sites
				reference_position += alignment.cigar.op_length(i);
				break;
			case BAM_CINS:
				score--;
				read_position += alignment.cigar.op_length(i);
				break;
			case BAM_CEQUAL:
				score += alignment.cigar.op_length(i);
				// fall through
			case BAM_CDIFF:
				reference_position += alignment.cigar.op_length(i);
				read_position += alignment.cigar.op_length(i);
				break;
			case BAM_CMATCH: // we need to check for mismatches ourselves
				for (unsigned int operation_i = 1; operation_i <= alignment.cigar.op_length(i); ++operation_i) {
					if (sequence[read_position] == assembly.at(alignment.contig)[reference_position])
						score++;
					reference_position++;
					read_position++;
				}
				break;
		}
	}
	return score;
}

int calculate_alignment_score(const mates_t& mates, const exon_annotation_index_t& exon_annotation_index, const assembly_t& assembly) {

	int score = calculate_segment_score(mates[MATE1], mates[MATE1].sequence, exon_annotation_index, assembly) +
	            calculate_segment_score(mates[MATE2], mates[MATE2].sequence, exon_annotation_index, assembly);

	if (mates.size() == 3) { // has a supplementary alignment
		score += calculate_segment_score(mates[SUPPLEMENTARY], (mates[SUPPLEMENTARY].strand == mates[SPLIT_READ].strand) ? mates[SPLIT_READ].sequence : dna_to_reverse_complement(mates[SPLIT_READ].sequence), exon_annotation_index, assembly);
		// penalize if the read is not split at a splice site
		if (!is_gap_at_splice_site((mates[SUPPLEMENTARY].strand == FORWARD) ? mates[SUPPLEMENTARY].end : mates[SUPPLEMENTARY].start, (mates[SUPPLEMENTARY].strand == FORWARD) ? DOWNSTREAM : UPSTREAM, mates[SUPPLEMENTARY].genes, exon_annotation_index) ||
		    !is_gap_at_splice_site((mates[SPLIT_READ].strand == FORWARD) ? mates[SPLIT_READ].start : mates[SPLIT_READ].end, (mates[SPLIT_READ].strand == FORWARD) ? UPSTREAM : DOWNSTREAM, mates[SPLIT_READ].genes, exon_annotation_index))
			score--;
	}

	return score;
}

// deterministically find the fusion with more supporting reads in a pair of fusions
bool fusion_has_more_support(const fusion_t* fusion, const fusion_t* current_best) {
	if (fusion == NULL) {
		return false;
	} else if (current_best == NULL) {
		return true;
	} else if (current_best->supporting_reads() != fusion->supporting_reads()) {
		return current_best->supporting_reads() < fusion->supporting_reads();
	} else if (fusion->gene1->is_protein_coding != current_best->gene1->is_protein_coding) { // preferrably pick protein-coding genes in case of ties
		return fusion->gene1->is_protein_coding;
	} else if (fusion->gene2->is_protein_coding != current_best->gene2->is_protein_coding) {
		return fusion->gene2->is_protein_coding;
	} else if (fusion->contig1 != current_best->contig1) { // all following rules are purely technical tie-breakers for deterministic behavior
		return fusion->contig1 < current_best->contig1;
	} else if (fusion->contig2 != current_best->contig2) {
		return fusion->contig2 < current_best->contig2;
	} else if (fusion->breakpoint1 != current_best->breakpoint1) {
		return fusion->breakpoint1 < current_best->breakpoint1;
	} else if (fusion->breakpoint2 != current_best->breakpoint2) {
		return fusion->breakpoint2 < current_best->breakpoint2;
	} else if (fusion->direction1 != current_best->direction1) {
		return fusion->direction1 < current_best->direction1;
	} else if (fusion->direction2 != current_best->direction2) {
		return fusion->direction2 < current_best->direction2;
	} else if (fusion->gene1->id != current_best->gene1->id) {
		return fusion->gene1->id < current_best->gene1->id;
	} else {
		return fusion->gene2->id < current_best->gene2->id;
	}
}


// get non-redundant supporting reads set
unordered_set<string> unique_supporting_reads_set(const fusion_t* fusion){
	unordered_set<string> unique_reads;
	for (auto chimeric_alignment = fusion->split_read1_list.begin(); chimeric_alignment != fusion->split_read1_list.end(); ++chimeric_alignment){
		if ((**chimeric_alignment).second.filter == FILTER_none || (**chimeric_alignment).second.filter == FILTER_multimappers)
			unique_reads.insert(strip_hi_tag_from_read_name((**chimeric_alignment).first));
	}
	for (auto chimeric_alignment = fusion->split_read2_list.begin(); chimeric_alignment != fusion->split_read2_list.end(); ++chimeric_alignment){
		if ((**chimeric_alignment).second.filter == FILTER_none || (**chimeric_alignment).second.filter == FILTER_multimappers)
			unique_reads.insert(strip_hi_tag_from_read_name((**chimeric_alignment).first));
	}
	for (auto chimeric_alignment = fusion->discordant_mate_list.begin(); chimeric_alignment != fusion->discordant_mate_list.end(); ++chimeric_alignment){
		if ((**chimeric_alignment).second.filter == FILTER_none || (**chimeric_alignment).second.filter == FILTER_multimappers)
			unique_reads.insert(strip_hi_tag_from_read_name((**chimeric_alignment).first));
	}
	return unique_reads;
}

// compared to "fusion_has_more_support", this function allows multiple alignment of the same score be kept for one read. 
// The first bool determine whether the alignment associated with the current fusion or the current best should be kept, the second bool determine whether the alignments of current fusion and current best should both be kept 
pair<bool, bool> fusion_has_more_non_redundant_support (const fusion_t* fusion, const fusion_t* current_best) {
	if (fusion == NULL) {
		return make_pair(false,false);
	} else if (current_best == NULL) {
		return make_pair(true,false);
	} else {
		unordered_set<string> current_best_set = unique_supporting_reads_set(current_best);
		unordered_set<string> fusion_set = unique_supporting_reads_set(fusion);
		if (fusion_set.size() != current_best_set.size() ){
			return make_pair(fusion_set.size() > current_best_set.size(),false);
		}else if (fusion->gene1->is_protein_coding != current_best->gene1->is_protein_coding) { // preferrably pick protein-coding genes in case of ties
			return make_pair(fusion->gene1->is_protein_coding,false);
		} else if (fusion->gene2->is_protein_coding != current_best->gene2->is_protein_coding) {
			return make_pair(fusion->gene2->is_protein_coding,false);
		} else {
			return make_pair(true,true);
		}
	}
}




// this function performs two tasks on multi-mapping reads:
// - when a read is assigned to multiple fusion candidates, is picks the candidate with the most supporting reads 
// - it selects the alignment with the highest alignment score
unsigned int filter_multimappers(chimeric_alignments_t& chimeric_alignments, fusions_t& fusions, const exon_annotation_index_t& exon_annotation_index, const assembly_t& assembly, disjoint_set& homolog_union) {

	// for each multi-mapper, find the fusion with the most supporting reads
	unordered_map<mates_t*,fusion_t*> most_supported_fusion;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		// cycle through split reads 1
		for (auto chimeric_alignment = fusion->second.split_read1_list.begin(); chimeric_alignment != fusion->second.split_read1_list.end(); ++chimeric_alignment) {
			fusion_t*& current_best = most_supported_fusion[&(**chimeric_alignment).second];
			if (fusion_has_more_support(&fusion->second, current_best))
				current_best = &fusion->second;
		}
		// cycle through split reads 2
		for (auto chimeric_alignment = fusion->second.split_read2_list.begin(); chimeric_alignment != fusion->second.split_read2_list.end(); ++chimeric_alignment) {
			fusion_t*& current_best = most_supported_fusion[&(**chimeric_alignment).second];
			if (fusion_has_more_support(&fusion->second, current_best))
				current_best = &fusion->second;
		}
		// cycle through discordant mates
		for (auto chimeric_alignment = fusion->second.discordant_mate_list.begin(); chimeric_alignment != fusion->second.discordant_mate_list.end(); ++chimeric_alignment) {
			fusion_t*& current_best = most_supported_fusion[&(**chimeric_alignment).second];
			if (fusion_has_more_support(&fusion->second, current_best))
				current_best = &fusion->second;
		}

		
	}



	// for each group of multi-mapping alignments, pick the one with the highest alignment score
	chimeric_alignments_t::iterator start_of_cluster = chimeric_alignments.begin();
	string read_name = (!chimeric_alignments.empty()) ? strip_hi_tag_from_read_name(start_of_cluster->first) : "";
	string next_read_name = read_name;
	string cluster_name = read_name;
	set<mates_t*> best_alignments;
	mates_t* best_alignment = NULL;
	int best_alignment_score = INT_MIN;
	
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); ; ++chimeric_alignment) {

		// the alignments are sorted by read number, so related alignments cluster together
		// check if we are still in the same cluster of multi-mapping reads
		read_name = next_read_name;
		if (cluster_name != read_name) {
			// this is the next cluster => go back over last cluster and discard all but the best multi-mapper
			if (!best_alignments.empty()){
				if(best_alignment==NULL){
					for (auto it = best_alignments.begin(); it != best_alignments.end(); ++it) {
						if(most_supported_fusion[*it]!=0){
							best_alignment = *it;
							break;
						}
    				}
					homolog_union.make_set(most_supported_fusion[best_alignment]);
				}
				for (chimeric_alignments_t::iterator chimeric_alignment2 = start_of_cluster; chimeric_alignment2 != chimeric_alignment; ++chimeric_alignment2){
					if (best_alignments.count(&chimeric_alignment2->second) > 0){
						if(most_supported_fusion[&chimeric_alignment2->second]!=0 && most_supported_fusion[best_alignment]!=0){
							homolog_union.make_set(most_supported_fusion[&chimeric_alignment2->second]);
							//std::cout << "after make_set: " << homolog_union.parent << '\n'
							//<< homolog_union.rank << std::endl;
							homolog_union.union_set(most_supported_fusion[&chimeric_alignment2->second],most_supported_fusion[best_alignment]);
							//std::cout << "after union_set: " << homolog_union.parent << '\n'
							//<< homolog_union.rank << std::endl;
						}
					}else{	
						if (chimeric_alignment2->second.filter == FILTER_none)
							chimeric_alignment2->second.filter = FILTER_multimappers;
					}
				}
			}


			// initialize variables for new cluster
			cluster_name = read_name;
			start_of_cluster = chimeric_alignment;
			best_alignments.clear();
			best_alignment = NULL;
			best_alignment_score = INT_MIN;
		}

		// abort once all alignments have been processed
		if (chimeric_alignment == chimeric_alignments.end())
			break;

		// peek into the next read
		// we can skip the calculation of alignment score, if it is a uniquely mapping read
		next_read_name = (next(chimeric_alignment) != chimeric_alignments.end()) ? strip_hi_tag_from_read_name(next(chimeric_alignment)->first) : "";
		if (start_of_cluster == chimeric_alignment && next_read_name != read_name)
			continue;

		// calculate alignment score and remember the best one
		int alignment_score = calculate_alignment_score(chimeric_alignment->second, exon_annotation_index, assembly);
		if (best_alignment_score < alignment_score) {
			best_alignments.clear();
			best_alignments.insert(&chimeric_alignment->second);
			best_alignment_score = alignment_score;
		} else if (best_alignment_score == alignment_score){ // when scores tie, pick the alignment that is associated with the fusion with more supporting reads
			if (fusion_has_more_non_redundant_support(most_supported_fusion[&chimeric_alignment->second], most_supported_fusion[best_alignment]).second){
				best_alignments.insert(&chimeric_alignment->second);
			}else if(fusion_has_more_non_redundant_support(most_supported_fusion[&chimeric_alignment->second], most_supported_fusion[best_alignment]).first){
				best_alignments.clear();
				best_alignments.insert(&chimeric_alignment->second);
			}
		}
			
	}

	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		//homologs[&fusion->second].clear();
		if (fusion->second.filter != FILTER_none || fusion->second.supporting_reads() == 0)
			continue;
		// scan split_read1_list
		for (auto chimeric_alignment = fusion->second.split_read1_list.begin(); chimeric_alignment != fusion->second.split_read1_list.end(); ++chimeric_alignment)
			if ((**chimeric_alignment).second.filter == FILTER_multimappers)
				if (fusion->second.split_reads1 > 0)
					fusion->second.split_reads1--;
		// scan split_read2_list
		for (auto chimeric_alignment = fusion->second.split_read2_list.begin(); chimeric_alignment != fusion->second.split_read2_list.end(); ++chimeric_alignment)
			if ((**chimeric_alignment).second.filter == FILTER_multimappers)
				if (fusion->second.split_reads2 > 0)
					fusion->second.split_reads2--;
		// scan discordant_mate_list
		for (auto chimeric_alignment = fusion->second.discordant_mate_list.begin(); chimeric_alignment != fusion->second.discordant_mate_list.end(); ++chimeric_alignment)
			if ((**chimeric_alignment).second.filter == FILTER_multimappers)
				if (fusion->second.discordant_mates > 0)
					fusion->second.discordant_mates--;

		if (fusion->second.supporting_reads() == 0){ // it was >0 before discarding multimapping reads
			fusion->second.filter = FILTER_multimappers; // now supporting_reads()==0 => the 'multimappers' filter discarded all supporting reads
		}else{
			if(fusion->second.gene1->name =="CSNK2A3" && fusion->second.gene2->name =="IGFBPL1"){
				
			}
		}
			
		

	}

	// count remaining fusions
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == FILTER_none)
			remaining++;
	return remaining;
}


