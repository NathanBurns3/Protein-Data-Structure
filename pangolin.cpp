#include "unit_test.h"
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <tuple>
#include <vector>
#include <string>


/* References:
 * https://bioboot.github.io/bimm143_W20/class-material/nw/
 * https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
 * https://www.ncbi.nlm.nih.gov/datasets/coronavirus/proteins/
 */

using direction = enum : uint8_t { diagonal = 0, horizontal = 1, vertical = 2 };
using NWAlignmentDifference = std::vector<std::vector<int>>;
using NWAlignmentDirection = std::vector<std::vector<direction>>;
using ProteinSequence = std::vector<char>;
using Comparison = std::vector<std::tuple<int, char, char>>;

/*
 * get_id_and_sequence
 *
 * It will get the identifier and list of amino acids from the
 * "next" protein (both as strings) in a FASTA-formatted file
 *
 * input: FASTA-formatted file, the id string, and sequence string
 * output: it will return false if there are no more identifier/sequence
 *         pairs to be read and will return true otherwise
 */
bool get_id_and_sequence(std::ifstream &file, std::string &id,
                         std::string &sequence) {
  if (!(!!std::getline(file, id))) {
    return false;
  }
  else if (!(!!std::getline(file, sequence))) {
    return false;
  }
  else {
    return true;
  }
}

/*
 * open_file
 *
 * It will open a file handle to the file named in the function parameter
 *
 * input: The filename as a string and an ifstream file variable
 * output: It will return false if it could not open the file and true otherwise
 */
bool open_file(std::string filename, std::ifstream &file) {
  file.open(filename);
  if (file.fail()) {
    return false;
  }
  else {
    return true;
  }
}

/*
 * string_to_protein_sequence
 *
 * It will convert a string to a sequence of nucleic acids
 *
 * input: a constant protein string variable
 * output: return the Protein Sequence data structure
 */
ProteinSequence string_to_protein_sequence(const std::string &protein_str) {
    ProteinSequence vect;
    for (int i = 0; i < protein_str.length(); i++) {
        vect.push_back(protein_str.at(i));
    }
    return vect;
}

/*
 * above_left
 *
 * takes a cell (row, col) and returns the value in the cell that is above/left of the inputed cell
 *
 * input: An NWAlignmentDifference variable named difference_vector
 *        and two int variables for the row and column
 * output: returns an int value for the selected cell
 */
int above_left(const NWAlignmentDifference &difference_vector, int row,
               int col) {
  /*
  If (row-1) and (col-1) are greater-than-or-equal-to zero, return the value in the cell at index (row-1), (col-1). Otherwise, if col equals 0, return -1 * row. Otherwise,
  return -1 * col. above_left will never be given a row/col pair that is less than 0 or greater than the bounds of difference_vector.
  */
  if ((row - 1) >= 0 && (col - 1) >= 0) {
      return difference_vector[row - 1][col - 1];
  }
  else if (col == 0) {
      return -1 * row;
  }
  else {
      return -1 * col;
  }
}

/*
 * above
 *
 * takes a cell (row, col) and returns the value in the cell that is above of the inputed cell
 *
 * input: An NWAlignmentDifference variable named difference_vector
 *        and two int variables for the row and column
 * output: returns an int value for the selected cell
 */
int above(const NWAlignmentDifference &difference_vector, int row, int col) {
  /*
   If (row-1) is greater-than-or-equal-to zero, return the value in the cell at index (row-1), col. Otherwise, return -1 * (col + 1).
   above will never be given a row/col pair that is less than 0 or greater than the bounds of difference_vector.
  */
  if ((row - 1) >= 0) {
      return difference_vector[row - 1][col];
  }
  else {
      return -1 * (col + 1);
  }
}

/*
 * left
 *
 * takes a cell (row, col) and returns the value in the cell that is left of the inputed cell
 *
 * input: An NWAlignmentDifference variable named difference_vector
 *        and two int variables for the row and column
 * output: returns an int value for the selected cell
 */
int left(const NWAlignmentDifference &difference_vector, int row, int col) {
  /*
  	If (col-1) is greater-than-or-equal-to zero, return the value in the cell at index row, (col-1) . Otherwise, return -1 * (row + 1).
    left will never be given a row/col pair that is less than 0 or greater than the bounds of difference_vector.
  */
  if ((col - 1) >= 0) {
      return difference_vector[row ][col - 1];
  }
  else {
      return -1 * (row + 1);
  }
}

void print_comparison(const Comparison &comparison,
                      const std::string &comparison_id, bool verbose = false) {

  std::cout << comparison_id << ": " << comparison.size();

  if (verbose) {
    for (auto i : comparison) {
      auto index = std::get<0>(i);
      auto first = std::get<1>(i);
      auto second = std::get<2>(i);
      std::cout << "(@" << index << ": " << first << ", " << second << "), ";
    }
  }
  std::cout << "\n";
}

template <typename T>
std::vector<T> reverse_vector(const std::vector<T> &to_reverse) {
  std::vector<T> forward;
  for (auto i : to_reverse) {
    forward.push_back(i);
  }
  return forward;
}

Comparison create_comparison(const ProteinSequence &top,
                             const ProteinSequence &bottom,
                             const NWAlignmentDirection &dirs) {
  int row = bottom.size() - 1;
  int col = top.size() - 1;
  Comparison full_result, diff_result;
  while (row >= 0 && col >= 0) {
    char first{0};
    char second{0};
    if (dirs[row][col] == diagonal) {
      first = bottom[row];
      second = top[col];
      row--;
      col--;
    } else if (dirs[row][col] == horizontal) {
      first = '-';
      second = top[col];
      col--;
    } else {
      first = bottom[row];
      second = '-';
      row--;
    }
    full_result.push_back(std::make_tuple(0, first, second));
  }
  full_result = reverse_vector(full_result);

  auto index = 0;
  for (auto i : full_result) {
    auto first = std::get<1>(i);
    auto second = std::get<2>(i);
    if (first != second) {
      diff_result.push_back(std::make_tuple(index, first, second));
    }
    index++;
  }
  return diff_result;
}

NWAlignmentDirection new_direction_vector(int rows, int cols) {
  std::vector<std::vector<direction>> new_direction_vector(
      rows, std::vector<direction>(cols));
  return new_direction_vector;
}

NWAlignmentDifference new_difference_vector(int rows, int cols) {
  NWAlignmentDifference new_diff_vector(rows, std::vector<int>(cols));
  return new_diff_vector;
}

int s(const char a, const char b) {
  if (a == 'X' || b == 'X') {
    return 0;
  }
  if (a == b) {
    return 1;
  }
  return -1;
}

NWAlignmentDirection needleman_wunsch(const ProteinSequence &top,
                                      const ProteinSequence &bottom) {
  auto bottom_length = bottom.size();
  auto top_length = top.size();
  auto directions = new_direction_vector(bottom_length, top_length);
  auto differences = new_difference_vector(bottom_length, top_length);
  for (int r = 0; r < bottom_length; r++) {
    for (int c = 0; c < top_length; c++) {
      int t, l, d, match, insert, del;
      t = above(differences, r, c);
      l = left(differences, r, c);
      d = above_left(differences, r, c);

      match = d + s(top[c], bottom[r]);
      del = t - 1;
      insert = l - 1;

      if (match > del && match > insert) {
        differences[r][c] = match;
        directions[r][c] = diagonal;
      } else if (del > insert) {
        differences[r][c] = del;
        directions[r][c] = vertical;
      } else {
        differences[r][c] = insert;
        directions[r][c] = horizontal;
      }
    }
  }
  return directions;
}

void run_unit_tests() {
  std::ifstream test_input_file;

  check_result(true, open_file("testing.fasta", test_input_file));
  check_result(true, test_input_file.is_open());

  std::string id, sequence;

  check_result(true, get_id_and_sequence(test_input_file, id, sequence));
  check_result(true, id == std::string{"id1"});
  check_result(true, sequence == std::string{"sequence1"});
  check_result(true, get_id_and_sequence(test_input_file, id, sequence));
  check_result(true, id == std::string{"id2"});
  check_result(true, sequence == std::string{"sequence2"});

  test_input_file.close();

  NWAlignmentDifference difference = new_difference_vector(2, 2);
  int counter{0};
  for (auto row : difference) {
    for (auto cell : row) {
      cell = counter++;
    }
  }

  check_result(0, above(difference, 1, 0));
  check_result(-2, above(difference, 0, 1));
  check_result(0, left(difference, 0, 1));
  check_result(-6, left(difference, 5, 0));
  check_result(0, above_left(difference, 1, 1));
  check_result(-4, above_left(difference, 0, 4));
  check_result(-3, above_left(difference, 3, 0));
}

int main() {
  //run_unit_tests();
  std::string id = "";
  std::string referenceSequence = "";
  std::string comparisonSequence = "";
  
  /*
  1. Open the file containing the reference protein sequence (reference.fasta)
    - check whether it opens (stop the program if the file cannot be opened)
  */
  std::ifstream referenceFile;
  if (!(open_file("reference.fasta", referenceFile)))
  {
    return 0;
  }
  
  /*
  2. Open the file containing the comparison protein sequences (comparison.fasta)
    - check whether it opens (stop the program if the file cannot be opened)
  */
  std::ifstream comparisonFile;
  if (!(open_file("comparison.fasta", comparisonFile)))
  {
    return 0;
  }

  /*
  3. Read the reference protein's id and sequence from reference.fasta using get_id_and_sequence
  */
  get_id_and_sequence(referenceFile, id, referenceSequence);

  /*
  4. Convert the reference sequence (which is a string) to a protein sequence (using string_to_protein_sequence)
  */
  ProteinSequence proteinReferenceSequence = string_to_protein_sequence(referenceSequence);

  /*
  5. As long as reading comparison id/sequence pairs from comparison.fasta (using get_id_and_sequence) succeeds,
    - convert the current comparison sequence (which is a string) to a protein sequence (using string_to_protein_sequence)
    - call the needleman_wunsch function to generate a NWAlignment Direction data structure
    - call create_comparison with the appropriate parameters (the results of the two previous functions) to generate a Comparison data structure
    - call print_comparison using the result of the previous step and the id of the current comparison sequence as parameters
  */
  while (get_id_and_sequence(comparisonFile, id, comparisonSequence)) {
    ProteinSequence proteinComparisonSequence = string_to_protein_sequence(comparisonSequence);
    NWAlignmentDirection direction = needleman_wunsch(proteinReferenceSequence, proteinComparisonSequence);
    Comparison comparison = create_comparison(proteinReferenceSequence, proteinComparisonSequence, direction);
    print_comparison(comparison, id);
  }
}