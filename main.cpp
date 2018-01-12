#include <iostream>

#include "polynom.h"


int main()
{

    polynom p;
    p.setFunc_name("f");
    p.setX("x");
    p+=4;
    p.add_parameter(0,1);
    p.add_parameter(1,1);
    p.add_parameter(2,1);
    p.add_parameter(3,1);
    std::cout << p << std::endl;
    std::cout << "Polynom size: " << p.degree() << std::endl;
    std::cout << polynom(p).compute_X_to_string(2) << std::endl;
    std::cout << p << std::endl;
    polynom d = p.derivative();
    std::cout << d << std::endl;
    polynom dtag = d.derivative();
    std::cout << dtag << std::endl;
    polynom dtag2 = dtag.derivative();
    std::cout << dtag2 << std::endl;
    polynom dtag3 = dtag2.derivative();
    std::cout << dtag3 << std::endl;

    std::cout << p << std::endl;
    std::cout << d.integral() << std::endl;

    polynom pkafold = p+d.integral();
    std::cout << pkafold.compute_X_to_string(2) << std::endl;

    std::cout << polynom(pkafold).compute_X_to_string(3) << std::endl;
    bool equality = polynom(pkafold)==pkafold;
    std::cout << equality << std::endl;

    std::cout << pkafold-polynom(pkafold) << std::endl;

    std::cout << (pkafold*d.integral()*p).compute_X_to_string(2) << std::endl;
    return 0;
}

