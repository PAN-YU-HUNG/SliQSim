#include "util_sim.h"
#include <string>
#include <fstream>
#include "Simulator.h"

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1(int length, int *reg)
{
    int one = 1, carry = 0;
    for (int i = 0; i < length; i++)
    {
        carry = reg[i] & one;
        reg[i] = reg[i] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1_start(int length, int *reg, int start)
{
    int one = 1, carry = 0;
    for (int i = start; i < length; i++)
    {
        carry = reg[i] & one;
        reg[i] = reg[i] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1_measure(int length, int *reg, int *order)
{
    int one = 1, carry = 0;
    for (int i = length - 1; i >= 0; i--)
    {
        carry = reg[order[i]] & one;
        reg[order[i]] = reg[order[i]] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int int_array_full_check(int length, int *reg)
{
    int check = 1;
    for (int i = 0; i < length; i++)
        check *= reg[i];

    return check;
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void replace_all(std::string& str, const std::string& search_str, const std::string& replace_str)
{
    std::string::size_type pos = 0;
    while ((pos = str.find(search_str, pos)) != std::string::npos)
    {
        str.replace(pos, search_str.size(), replace_str);
        pos += replace_str.size();
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
std::string find_in_yaml(std::istream& yaml, const std::string& key)
{
    // std::cout<<"key: "<<key<<std::endl;
    std::string res;
    std::string line;
    while(std::getline(yaml, line)){
        // std::cout<<"line: "<<line<<std::endl;
        auto line_key = line.substr(0, line.find(":"));
        // std::cout<<"line_key: "<<line_key<<std::endl;
        if (line_key == key){
             res = line.substr(line.find(":") + 1);
             break;
        }
    }
    if(res[0] == ' '){
        res = res.substr(1);
    }
    return res;
}

/**Function*************************************************************

  Synopsis    []

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
double parse_theta(std::string theta){
  try{
    auto pi_pos = theta.find("pi");
    bool with_pi = pi_pos != std::string::npos;
    double scaler = with_pi ? PI: 1;
    if(with_pi){
        theta.replace(pi_pos, 2, "");
    }
    auto slash_pos = theta.find("/");
    if(slash_pos != std::string::npos){
      std::string numerator, denominator;
      numerator = theta.substr(0, slash_pos);
      if(numerator.empty()){
        numerator = "1";
      } else if(numerator == "-"){
        numerator = "-1";
      }
      denominator = theta.substr(slash_pos + 1);
      return std::stod(numerator) / std::stod(denominator) * scaler;
    }
    return std::stod(theta) * scaler;
    } catch(std::invalid_argument& e){
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << "Invalid theta: " << theta << std::endl;
        throw e;
    }
}

/**Function*************************************************************

  Synopsis    [normalize to [-pi, pi]]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
double normalize_theta(double theta){
    while(theta > PI){
        theta -= 2 * PI;
    }
    while(theta < -PI){
        theta += 2 * PI;
    }
    return theta;
}