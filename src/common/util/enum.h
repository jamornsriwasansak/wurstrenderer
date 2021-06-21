#pragma once

#define DECLARE_ENUM_OPERATORS(Type) \
Type operator|(const Type a, const Type b)\
{ return static_cast<Type>(static_cast<unsigned char>(a) | static_cast<unsigned char>(b)); }\
Type operator&(const Type a, const Type b)\
{ return static_cast<Type>(static_cast<unsigned char>(a) & static_cast<unsigned char>(b)); }\
bool operator==(const Type a, const Type b)\
{ return static_cast<unsigned char>(a) == static_cast<unsigned char>(b); }\
bool operator!=(const Type a, const Type b)\
{ return !(a == b); }\
 bool operator==(const Type a, unsigned char b)\
{ return static_cast<unsigned char>(a) == b; }\
bool operator!=(const Type a, unsigned char b)\
{ return !(a == b); }
