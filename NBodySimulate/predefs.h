#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdint>
#include <chrono>

struct Vector2 {
    double x = 0, y = 0;

    Vector2()
        : x(0)
        , y(0) {}

    Vector2(double x, double y)
        : x(x)
        , y(y) {}

    Vector2 operator+(const Vector2& other) const noexcept {
        return { x + other.x, y + other.y };
    }

    Vector2 operator-(const Vector2& other) const noexcept {
        return { x - other.x, y - other.y };
    }

    Vector2 operator*(double num) const noexcept {
        return { x * num, y * num };
    }

    Vector2 operator/(double num) const {
        return { x / num, y / num };
    }

    Vector2& operator+=(const Vector2& other) noexcept {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2& operator-=(const Vector2& other) noexcept {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2& operator*=(const double num) noexcept {
        x *= num;
        y *= num;
        return *this;
    }

    Vector2& operator/=(const double num) {
        x /= num;
        y /= num;
        return *this;
    }

    bool operator==(const Vector2& other) const noexcept {
        return x == other.x && y == other.y;
    }

    bool operator!=(const Vector2& other) const noexcept {
        return x != other.x || y != other.y;
    }

    bool at_rd_of(const Vector2& other) const noexcept {
        return x > other.x && y > other.y;
    }

    double len() const {
        return std::pow(x * x + y * y, 0.5);
    }

    double len_square() const {
        return x * x + y * y;
    }

    Vector2 unit() const {
        return *this / len();
    }
};

inline Vector2 operator*(const double num, const Vector2& v) noexcept {
    return v * num;
}

inline std::ostream& operator<<(std::ostream& out, const Vector2& vec) {
    out << "(" << vec.x << ", " << vec.y << ")";
    return out;
}

struct Body {
    Vector2 pos, v;
    double mass;

    Body()
        : pos(.0, .0)
        , v(.0, .0)
        , mass(.0) {}

    Body(const Vector2& pos, double mass)
        : pos(pos)
        , mass(mass) {
    }

    Body(const Vector2& pos, const Vector2& velocity, double mass)
        : pos(pos)
        , v(velocity)
        , mass(mass) {
    }

    bool operator==(const Body& other) const noexcept {
        return mass == other.mass && pos == other.pos;
    }

    bool operator!=(const Body& other) const noexcept {
        return !(*this == other);
    }

	void move_under_force(Vector2 force, double seconds) {
		uniform_linear_move_for(seconds);
		add_force_to_speed(force, seconds);
	}

    void add_force_to_speed(Vector2 force, double seconds) {
        // v = v0 + at = v0 + (F/m)t
        // v += F(t/m)

        this->v += force * (seconds / mass);
    }

	void uniform_linear_move_for(double seconds) {
		pos += this->v * seconds;
	}
};

inline std::ostream& operator<<(std::ostream& out, const Body& body) {
    out << "Body { " << body.pos << ", " << body.mass << "kg }";
    return out;
}

struct BHNode {
    bool internal = false;
    int16_t parent = -1;

    // Virtual or Actual Body
    //（当 internal == true 时，为 virtual，否则为 actual）
    Body vaBody;
    int16_t child[4] = { -1, -1, -1, -1 };
};

struct Settings {
    size_t bodyNum = SIZE_MAX;
};

enum class Quadrant {
    RU = 0,
    LU,
    LD,
    RD
};
