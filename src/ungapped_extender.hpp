/*
 * ungapped_extender.hpp
 *
 *  Created on: 2012/11/12
 *      Author: shu
 */

#ifndef UNGAPPED_EXTENDER_HPP
#define UNGAPPED_EXTENDER_HPP

#include "alphabet_coder.hpp"
#include "edit_blocks.hpp"
#include "score_matrix.hpp"
#include <stdint.h>

class UngappedExtender {
public:
	UngappedExtender(const ScoreMatrix &score_matrix_, const AlphabetCoder::Code delimiter_)
			:number_letters(score_matrix_.GetNumberLetters()),sm(score_matrix_.GetMatrix()),delimiter(delimiter_){};

	int UngappedExtend(const AlphabetCoder::Code *query,
										const uint32_t query_position,
										const AlphabetCoder::Code *database,
										const uint32_t database_position,
										const int cutoff, const int trigger) const;

	int ExtendOneSideScoreOnly1(const AlphabetCoder::Code *sequence0,
			const AlphabetCoder::Code *sequence1, int cutoff) const;
    //__attribute__((noinline));
	int ExtendOneSideScoreOnly2(const AlphabetCoder::Code *sequence0,
			const AlphabetCoder::Code *sequence1, int cutoff) const;

	/*static bool ExtendOneSide(const AlphabetCoder::Code *sequence0,
			const AlphabetCoder::Code *sequence1,
			const AlphabetCoder::Code sequence_delimiter, bool reversed,
			const ScoreMatrix &score_matrix, int cutoff, int *best_score,
			int *best_sequence0_position, int *best_sequence1_position,
      EditBlocks *edit_blocks); //__attribute__((noinline));
*/
	//UngappedExtender();
	//virtual ~UngappedExtender();
private:
	const uint32_t number_letters;
	const int *sm;
	AlphabetCoder::Code delimiter;
};

inline int UngappedExtender::UngappedExtend(const AlphabetCoder::Code *query,
									const uint32_t query_position,
									const AlphabetCoder::Code *database,
									const uint32_t database_position,
									const int cutoff, const int trigger) const{
	
	int sum_score = ExtendOneSideScoreOnly1(query + query_position,
											database + database_position, cutoff);

	if(sum_score > trigger) return sum_score;

	sum_score += ExtendOneSideScoreOnly2(query + query_position -1,
								database + database_position -1,cutoff);

	return sum_score;
}

inline int UngappedExtender::ExtendOneSideScoreOnly1(
		const AlphabetCoder::Code *sequence0,
		const AlphabetCoder::Code *sequence1, int cutoff) const {

	int threshold = -cutoff;
	int position = 0;
	int b_score = 0;
	int score = 0;

	do {
		//cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score <<endl;
		if (sequence1[position] == delimiter
				|| sequence0[position] == delimiter) {
			break;
		}
		score += sm[sequence1[position] * number_letters
				+ sequence0[position]];
		++position;
		if (score > 0) {
			do {
				b_score += score;
				//cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score <<endl;
				if (sequence1[position] == delimiter
						|| sequence0[position] == delimiter) {
					break;
				}
				score = sm[sequence1[position] * number_letters
						+ sequence0[position]];
				++position;
			} while (score > 0);
		}
	} while (score > threshold);

	return b_score;
}

inline int UngappedExtender::ExtendOneSideScoreOnly2(
		const AlphabetCoder::Code *sequence0,
		const AlphabetCoder::Code *sequence1, int cutoff) const {

	int threshold = -cutoff;
	int position = 0;
	int b_score = 0;
	int score = 0;

	do {
		//cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score <<endl;
		if (sequence1[position] == delimiter
				|| sequence0[position] == delimiter) {
			break;
		}
		score += sm[sequence1[position] * number_letters
				+ sequence0[position]];
		--position;
		if (score > 0) {
			do {
				b_score += score;
				//cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score <<endl;
				if (sequence1[position] == delimiter
						|| sequence0[position] == delimiter) {
					break;
				}
				score = sm[sequence1[position] * number_letters
						+ sequence0[position]];
				--position;
			} while (score > 0);
		}
	} while (score > threshold);

	return b_score;
}
/*
inline bool UngappedExtender::ExtendOneSide(
		const AlphabetCoder::Code *sequence0,
		const AlphabetCoder::Code *sequence1,
		const AlphabetCoder::Code sequence_delimiter, bool reversed,
		const ScoreMatrix &score_matrix, int cutoff, int *best_score,
		int *best_sequence0_position, int *best_sequence1_position,
		EditBlocks *edit_blocks) {
	uint32_t number_letters = score_matrix.GetNumberLetters();
	const int *sm = score_matrix.GetMatrix();
	int increment = 0;
	if (reversed) {
		increment = -1;
	} else {
		increment = 1;
	}
	*best_sequence0_position = -increment;
	*best_sequence1_position = -increment;
	int threshold = -cutoff;
	int position = 0;
	int b_position = 0;
	int b_score = 0;
	if (sequence1[position] != sequence_delimiter
			&& sequence0[position] != sequence_delimiter) {
		int score = 0;
		do {
			//cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score <<endl;
			if (sequence1[position] == sequence_delimiter
					|| sequence0[position] == sequence_delimiter) {
				break;
			}
			score += sm[sequence1[position] * number_letters
					+ sequence0[position]];
			position += increment;
			if (score > 0) {
				do {
					b_score += score;
					b_position = position;
					//cout << position << " " <<(int)sequence0[position] << " " << (int)sequence1[position] << " "  <<  score <<endl;
					if (sequence1[position] == sequence_delimiter
							|| sequence0[position] == sequence_delimiter) {
						break;
					}
					score = sm[sequence1[position] * number_letters
							+ sequence0[position]];
					position += increment;
				} while (score > 0);
			}
		} while (score > threshold);
	}
	b_position = b_position - increment;
	*best_sequence0_position = b_position;
	*best_sequence1_position = b_position;
	*best_score = b_score;

	if (edit_blocks != NULL) {
		int edit_length = (b_position + increment) * increment;
		edit_blocks->Add(EditBlocks::kSubstitution, edit_length);
	}

	return true;
}*/

#endif /* UNGAPPED_EXTENDER_HPP */
