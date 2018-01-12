
#include "polynom.h"
#include <sstream>
#include <algorithm>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
polynom::polynom()
{
    this->p = new std::vector<double>();
    this->tags = 0;
    this->func_name = "f";
    this->x = "x";
    this->intgral_signed=false;
}

polynom::polynom(const polynom &other)
{
    this->p = new std::vector<double>(other.p->size());
    #pragma omp parallel for
    for(int i=0; i<other.p->size(); i++)
    {
        double val = other[i];
        (*this->p)[i] = val;
    }
    this->tags = 0;
    this->func_name = std::string(other.func_name);
    this->x = std::string(other.x);
    this->intgral_signed = other.intgral_signed;
}

polynom::~polynom()
{
    delete this->p;
    this->func_name.clear();
    this->x.clear();
    this->intgral_signed = false;
}

std::string polynom::toString() const
{
    std::string return_value;
    if(this->intgral_signed)
    {
        return_value+="[";
        return_value+=u8"∫";
        return_value+="] ";
    }
    return_value+=this->func_name;
    #pragma omp parallel for
    for(auto i=0; i<this->tags;i++)
        return_value+="'";
    return_value+=+"("+this->x+") = ";
    auto size = return_value.size();
    int index=0;
    bool first_flag = false;
    for(double value : *this->p)
    {
        if(value==0)
        {
            index++;
            continue;
        }
        if(first_flag)
            return_value+= " + ";
        std::string value_string;
        if(value == (int) value)
            value_string = std::to_string((int) value);
        else
            value_string = std::to_string(value);
        if(value != 1)
        {
             if(index==0)
             {
                 return_value += value_string ;
             }
             else
             {
                 return_value += value_string + "*" ;
             }
        }
        else
        {
            if(index==0)
                return_value += value_string ;
        }
        if(index!=0)
        {
            return_value += this->x;
            if(index!=1)
            {
                return_value += "^" + std::to_string(index);
            }
        }
        index++;
        first_flag=true;
    }
    if(size == return_value.size())
        return_value+="0";
    if(this->intgral_signed)
        return_value+=" + c";
    return return_value;
}

std::ostream &operator<<(std::ostream &os, const polynom &polynom1)
{

    os << polynom1.toString();
    return os;
}


bool polynom::operator==(const polynom &rhs) const
{
    if(this->p->size() != rhs.p->size())
        return false;
    bool return_value = true;
    #pragma omp parallel for
    for(auto i= static_cast<int>(this->p->size() - 1); i >= 0; i--)
    {
        if(this->p->at(static_cast<unsigned long>(i)) != p->at(static_cast<unsigned long>(i)))
            return_value = false;
    }
    return return_value;
}

bool polynom::operator!=(const polynom &rhs) const {
    return !(rhs == *this);
}

bool polynom::operator<(const polynom &rhs) const {
    if(this->p->size() < rhs.p->size())
        return true;
    bool return_value = false;
    #pragma omp parallel for
    for(auto i= static_cast<int>(this->p->size() - 1); i >= 0; i--)
    {
        if(this->p->at(static_cast<unsigned long>(i)) < p->at(static_cast<unsigned long>(i)))
            return_value = true;
    }
    return return_value;
}

bool polynom::operator>(const polynom &rhs) const {
    return !(this->operator==(rhs) || this->operator<(rhs));
}

bool polynom::operator<=(const polynom &rhs) const {
    return !(rhs < *this);
}

bool polynom::operator>=(const polynom &rhs) const {
    return !(*this < rhs);
}

double polynom::compute_X(double x) const
{
    double sum =0;
    double arr[this->p->size()];
    #pragma omp parallel for
    for(int i=0;i<this->p->size();i++)
        arr[i] = 0;
    #pragma omp parallel for
    for(int i=0; i<this->p->size();i++)
    {
        double value = this->p->at(static_cast<unsigned long>(i));
        if(value != 0)
            arr[i]=value*std::pow(x,i);
    }

    for(int i=0; i<this->p->size();i++)
        sum+=arr[i];
    return sum;
}

polynom &polynom::operator++()
{
    if(this->p->empty())
        this->p->push_back(1);
    else
        this->p->assign(0,this->p->at(0)+1);
    return *this;
}

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
polynom polynom::operator++(int op)
#pragma clang diagnostic pop
{
    polynom& copy(*this);
    ++(*this);
    return copy;
}

polynom& polynom::operator+=(int op)
{
    if(this->p->empty())
        this->p->push_back(op);
    else
        this->p->assign(0,this->p->at(0)+op);
    return *this;
}

polynom& polynom::operator-=(int op) {
    return this->operator+=(-op);
}

void polynom::add_parameter(int index, double parameter)
{
    if(index < 0)
        return;
    auto maximum = (unsigned long) std::max(index+1,(int) p->size());
    this->p->resize(maximum);
    (*this->p)[static_cast<unsigned long>(index)] = parameter;
}

int polynom::degree() const
{
    auto size = static_cast<int>(this->p->size() - 1);
    while((*this->p)[static_cast<unsigned long>(size)]==0 && size>=0)
        size--;
    return size;
}

polynom polynom::derivative() const
{
    polynom d;
    d.p->resize(this->p->size());
    #pragma omp parallel for
    for(int i=1; i<this->p->size(); i++)
    {
        double val = i*((*this)[i]);
        (*d.p)[i-1] = val;
    }
    d.setFunc_name(this->func_name);
    d.setX(this->x);
    d.tags = this->tags+1;
    d.intgral_signed = false;
    return d;
}

polynom polynom::integral() const
{
    polynom d;
    d.p->resize(this->p->size());
    #pragma omp parallel for
    for(auto i=0; i<this->p->size();i++)
    {
        double val = (*this)[i];
        (*d.p)[i] = val/(i+1);
    }
    d.p->insert( d.p->begin(), d.p->back() );
    d.add_parameter(0,0);
    d.setFunc_name(this->func_name);
    d.setX(this->x);
    d.tags = 0;
    d.intgral_signed = true;
    return d;
}

std::string polynom::compute_X_to_string(double x) const
{
    std::string return_value;
    if(this->intgral_signed)
    {
        return_value+="[";
        return_value+=u8"∫";
        return_value+="] ";
    }
    return_value+=this->func_name;
    #pragma omp parallel for
    for(auto i=0; i<this->tags;i++)
        return_value+="'";
    if(x == (int) x)
        return_value+="(" + std::to_string((int) x) + ") = ";
    else
        return_value+= "("+ std::to_string(x)+") = ";
    auto val = this->compute_X(x);
    if(val == (int) val)
        return_value += std::to_string((int) val);
    else
        return_value += std::to_string(val);
    if(this->intgral_signed)
        return_value+=" + c";
    return return_value; // NOLINT
}

const std::string &polynom::getFunc_name() const {
    return func_name;
}

void polynom::setFunc_name(const std::string &func_name) {
    polynom::func_name = std::string(func_name);
}

const std::string &polynom::getX() const {
    return x;
}

void polynom::setX(const std::string &x) {
    polynom::x = std::string(x);
}

int polynom::operator[](int x) const
{
    if(x>=this->p->size() || x<0)
        return 0;
    else
        return static_cast<int>((*this->p)[static_cast<unsigned long>(x)]);
}

polynom polynom::operator+(const polynom &other)
{
    polynom three;
    auto maximum = (unsigned long) std::max((int) other.p->size(),(int) p->size());
    three.p->resize(maximum);
    #pragma omp parallel for
    for(int n=0; n<maximum; ++n)
    {
        (*three.p)[n] = (*this)[n] + other[n];
    }
    three.setFunc_name(this->func_name);
    three.setX(this->x);
    three.intgral_signed=false;
    three.tags=0;
    return three;
}

polynom polynom::operator-(const polynom &other)
{
    polynom three;
    auto maximum = (unsigned long) std::max((int) other.p->size(),(int) p->size());
    three.p->resize(maximum);
    #pragma omp parallel for
    for(int n=0; n<maximum; ++n)
    {
        (*three.p)[n] = (*this)[n] - other[n];
    }
    three.setFunc_name(this->func_name);
    three.setX(this->x);
    three.intgral_signed=false;
    three.tags=0;
    return three;
}

polynom polynom::operator*(const polynom &other)
{
    polynom three;
    for (int i=0;i<this->p->size();i++)
    {
        for (int j=0;j<other.p->size();j++)
        {
            three.add_parameter(i+j,three[i+j] + ((*this)[i] * other[j]));
        }
    }
    three.setFunc_name(this->func_name);
    three.setX(this->x);
    three.intgral_signed=false;
    three.tags=0;
    return three;
}

polynom polynom::operator+=(const polynom &other)
{
    return (*this)+other;
}

polynom polynom::operator-=(const polynom &other)
{
    return (*this)-other;
}

polynom polynom::operator*=(const polynom &other) {
    return (*this)*other;
}




#pragma clang diagnostic pop