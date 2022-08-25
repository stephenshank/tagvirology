#include <string>
#include <filesystem>
 
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
 
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
 
    for (auto & [seq, id, qual] : file_in)
    {
        seqan3::debug_stream << "ID:     " << id << '\n';
        seqan3::debug_stream << "SEQ:    " << seq << '\n';
        seqan3::debug_stream << "Empty Qual." << qual << '\n'; // qual is empty for FASTA files
    }
 
    return 0;
}
