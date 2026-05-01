#include "cas/bigint.hpp"
#include <iostream>
int main() {
    cas::BigInt b(-4);
    std::cout << "decimal: " << b.decimal() << std::endl;
    return 0;
}