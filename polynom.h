
#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
#ifndef PREPARE_FOR_EXAM_POLYNOM_H
#define PREPARE_FOR_EXAM_POLYNOM_H

#include <vector>
#include <iostream>
#include <cmath>
#include <omp.h>

class polynom {
private:
    std::vector<double >* p;
    int tags = 0;
    bool intgral_signed;
    std::string func_name;
    std::string x;
public:
    polynom();
    polynom(const polynom& other);

    friend std::ostream &operator<<(std::ostream &os, const polynom &polynom1);

    bool operator==(const polynom &rhs) const;

    bool operator!=(const polynom &rhs) const;

    bool operator<(const polynom &rhs) const;

    bool operator>(const polynom &rhs) const;

    bool operator<=(const polynom &rhs) const;

    bool operator>=(const polynom &rhs) const;

    int operator[](int x) const;

    const std::string &getFunc_name() const;

    void setFunc_name(const std::string &func_name);

    const std::string &getX() const;

    void setX(const std::string &x);

    polynom& operator++();
    polynom operator++(int op);

    polynom& operator+=(int op);
    polynom& operator-=(int op);

    polynom operator+=(const polynom& other);
    polynom operator-=(const polynom& other);

    polynom operator+(const polynom& other);
    polynom operator-(const polynom& other);

    polynom operator*(const polynom& other);
    polynom operator*=(const polynom& other);


    void add_parameter(int index, double parameter);

    int degree() const;

    double compute_X(double x) const;
    std::string compute_X_to_string(double x) const;

    polynom derivative() const;

    polynom integral() const ;

    std::string toString() const;

    virtual ~polynom();
};


#endif //PREPARE_FOR_EXAM_POLYNOM_H

#pragma clang diagnostic pop