/******************************************************
* Constant/Variable polymorphing helper
*
* (c) 2021 Alexander Petrov-Savchenko
* e-mail: axp@soft-amp.com
*
******************************************************/
#pragma once 
#ifndef __DSY_CONSTVAR_H__
#define __DSY_CONSTVAR_H__

namespace daisysp
{

template<typename T, T init = 0, bool fixed = false>
class ConstVar
{
public:
    ConstVar()
    :   value_(init)
    {
    }

    ConstVar(const T& val)
    {
        value_ = val;
    }
    // void Set(const T& value)
    // {
    //     value_ = value;
    // }
    operator T&() noexcept 
    {
        return value_;
    }
    operator T() const noexcept
    {
        return value_;
    }
protected:
    T value_;
};


template<typename T, T value> 
struct ConstVar<T, value, true>
{
public:
    ConstVar()
    {
    }
  //  static void Set(T __unused new_value)
  //  {        
   // }
    static constexpr T Get()
    {
        return value;
    }
};

}  // namespace daisysp
#endif // __CONSTVAR_H__