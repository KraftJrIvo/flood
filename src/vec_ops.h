#pragma once

#include <iostream>
#include "raylib.h"
#include "raymath.h"

inline std::ostream& operator<<(std::ostream& os, const Vector2& v) {
    os << v.x << ' ' <<  v.y;
    return os;
}

inline bool operator== (const Vector2& v1, const Vector2& v2) {
    return v1.x == v2.x && v1.y == v2.y;
}

inline Vector2 operator+ (const Vector2& v1, const Vector2& v2) {
    return { v1.x + v2.x, v1.y + v2.y };
}

inline Vector2 operator- (const Vector2& v1, const Vector2& v2) {
    return { v1.x - v2.x, v1.y - v2.y };
}

inline Vector2 operator- (const Vector2& v) {
    return { -v.x, -v.y };
}

inline void operator+=(Vector2& v1, const Vector2& v2) {
    v1.x += v2.x;
    v1.y += v2.y;
}

inline void operator-=(Vector2& v1, const Vector2& v2) {
    v1.x -= v2.x;
    v1.y -= v2.y;
}

inline Vector2 operator* (const Vector2& v, const float& coeff) {
    return { coeff * v.x, coeff * v.y };
}

inline Vector2 operator/ (const Vector2& v, const float& coeff) {
    return { v.x / coeff, v.y / coeff };
}

inline Vector2 operator* (const float& coeff, const Vector2& v) {
    return v * coeff;
}

inline Vector2 operator/ (const float& coeff, const Vector2& v) {
    return v / coeff;
}

inline Vector2 operator* (const Vector2& v1, const Vector2& v2) {
    return { v1.x * v2.x, v1.y * v2.y };
}

inline std::ostream& operator<<(std::ostream& os, const Vector3& v) {
    os << v.x << ' ' <<  v.y << ' ' <<  v.z;
    return os;
}

inline bool operator== (const Vector3& v1, const Vector3& v2) {
    return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

inline Vector3 operator+ (const Vector3& v1, const Vector3& v2) {
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

inline Vector3 operator- (const Vector3& v1, const Vector3& v2) {
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

inline Vector3 operator- (const Vector3& v) {
    return { -v.x, -v.y, -v.z };
}

inline void operator+=(Vector3& v1, const Vector3& v2) {
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
}

inline void operator-=(Vector3& v1, const Vector3& v2) {
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
}

inline Vector3 operator* (const Vector3& v, const float& coeff) {
    return { coeff * v.x, coeff * v.y, coeff * v.z };
}

inline Vector3 operator/ (const Vector3& v, const float& coeff) {
    return { v.x / coeff, v.y / coeff, v.z / coeff };
}

inline Vector3 operator* (const float& coeff, const Vector3& v) {
    return v * coeff;
}

inline Vector3 operator/ (const float& coeff, const Vector3& v) {
    return v / coeff;
}

inline Vector3 operator* (const Vector3& v1, const Vector3& v2) {
    return { v1.x * v2.x, v1.y * v2.y, v1.z * v2.z };
}

inline float& at(Vector3& v, int idx) {
    return (idx == 0) ? v.x : (idx == 1) ? v.y : v.z;
}

inline float& at(Matrix& m, int idx1, int idx2) {
    return (idx1 == 0) ? ((idx2 == 0) ? m.m0 : (idx2 == 1) ? m.m4 : (idx2 == 2) ? m.m8 : (idx2 == 3) ? m.m12 : m.m0) :
           (idx1 == 1) ? ((idx2 == 0) ? m.m1 : (idx2 == 1) ? m.m5 : (idx2 == 2) ? m.m9 : (idx2 == 3) ? m.m13 : m.m0) :
           (idx1 == 2) ? ((idx2 == 0) ? m.m2 : (idx2 == 1) ? m.m6 : (idx2 == 2) ? m.m10 : (idx2 == 3) ? m.m14 : m.m0) :
                         ((idx2 == 0) ? m.m3 : (idx2 == 1) ? m.m7 : (idx2 == 2) ? m.m11 : (idx2 == 3) ? m.m15 : m.m0);
}

inline float at(const Matrix& m, int idx1, int idx2) {
    return (idx1 == 0) ? ((idx2 == 0) ? m.m0 : (idx2 == 1) ? m.m4 : (idx2 == 2) ? m.m8 : (idx2 == 3) ? m.m12 : m.m0) :
        (idx1 == 1) ? ((idx2 == 0) ? m.m1 : (idx2 == 1) ? m.m5 : (idx2 == 2) ? m.m9 : (idx2 == 3) ? m.m13 : m.m0) :
        (idx1 == 2) ? ((idx2 == 0) ? m.m2 : (idx2 == 1) ? m.m6 : (idx2 == 2) ? m.m10 : (idx2 == 3) ? m.m14 : m.m0) :
        ((idx2 == 0) ? m.m3 : (idx2 == 1) ? m.m7 : (idx2 == 2) ? m.m11 : (idx2 == 3) ? m.m15 : m.m0);
}