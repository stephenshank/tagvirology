#include <string>
#include <filesystem>
 
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
 
int main(int argc, char * argv[])
{
    std::string input{};
    std::string prefix{};
    seqan3::argument_parser parser("barcode", argc, argv);
    parser.add_positional_option(input, "Path to file");
    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
        return 0;
    }

    bool is_absolute_path = input[0] == '/';
    if (is_absolute_path) {
      prefix = "";
    } else {
      prefix = std::filesystem::current_path().string() + "/";
    }
    seqan3::sequence_file_input file_in{prefix + input};

    using namespace seqan3::literals;
    seqan3::dna4_vector constant_region_1 = "AGTACGTACGAGTC"_dna4;
    seqan3::dna4_vector constant_region_2 = "GTACTCGCAGTAGTC"_dna4;

    auto config = seqan3::align_cfg::method_global{
      seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
      seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
      seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
      seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}
    } |
      seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}};

    int i;
    seqan3::debug_stream << "id\talignment_1\tstart_1\tend_1\tbarcode_1\t"
      << "alignment_2\tstart_2\tend_2\tbarcode_2\n";
    for (auto & [seq, id, qual] : file_in)
    {
        seqan3::debug_stream << id << '\t';

        auto results_1 = seqan3::align_pairwise(std::tie(seq, constant_region_1), config);
        auto & res1 = *results_1.begin();
        seqan3::debug_stream << res1.alignment() << '\t'
          << res1.sequence1_begin_position()
          << "\t" <<
          res1.sequence1_end_position() << "\t";
        for(i=8; i>0; i--) {
          seqan3::debug_stream << seq[res1.sequence1_begin_position() - i - 1];
        }
        seqan3::debug_stream << "\t";

        auto results_2 = seqan3::align_pairwise(std::tie(seq, constant_region_2), config);
        auto & res2 = *results_2.begin();
        seqan3::debug_stream << res2.alignment() << '\t'
          << res2.sequence1_begin_position()
          << "\t" <<
          res2.sequence1_end_position() << "\t";
        for(i=8; i>0; i--) {
          seqan3::debug_stream << seq[res2.sequence1_begin_position() - i - 1];
        }
        seqan3::debug_stream << "\t\n";
    }
 
    return 0;
}
