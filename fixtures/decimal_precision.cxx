#include <correction.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

int main(int argc, char** argv) {
  if (argc != 3) return 2;
  std::ifstream input(argv[1]);
  std::ostringstream text;
  text << input.rdbuf();
  const auto from_file = correction::CorrectionSet::from_file(argv[1]);
  const auto from_gzip = correction::CorrectionSet::from_file(argv[2]);
  const auto from_string = correction::CorrectionSet::from_string(text.str().c_str());
  const double expected = 0x1.b6e246f39336ep+10;
  const double file_value = from_file->at("decimal")->evaluate({});
  const double string_value = from_string->at("decimal")->evaluate({});
  const double gzip_value = from_gzip->at("decimal")->evaluate({});
  std::cout << std::hexfloat << "expected " << expected << '\n'
            << "file " << file_value << '\n'
            << "gzip " << gzip_value << '\n'
            << "string " << string_value << '\n';
  return file_value == expected && string_value == expected && gzip_value == expected ? 0 : 1;
}
