#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "raylib.h"
#include "raymath.h"

#include "vec_ops.h"

#define INPUT_IMG "res/test3.png"
#define GRAVITY Vector2{0.0f, 0.1f}
#define PRESSURE 0.222f
#define EPS 1e-3f

#define SCALE 4.0f

struct FloodCell {
    float mass;
    Vector2 vel;
    Vector3 color;

    void operator+=(const FloodCell& b) {
        vel = vel * mass / (mass + b.mass) + b.vel * b.mass / (mass + b.mass);
        color = color * mass / (mass + b.mass) + b.color * b.mass / (mass + b.mass);
        mass += b.mass;
    }
    FloodCell operator*(float coeff) const {
        FloodCell fc;
        fc.vel = vel;// * coeff;
        fc.mass = mass * coeff;
        fc.color = color;// * coeff;
        return fc;
    }
};

typedef std::vector<std::vector<FloodCell>> FloodField;

FloodField makeField(std::string imagePath) {
    auto image = LoadImage(imagePath.c_str());
    FloodField field(image.height, std::vector<FloodCell>(image.width, FloodCell{0.0f, Vector2Zero(), Vector3Zero()}));
    for (int i = 0; i < image.height; ++i) {
        for (int j = 0; j < image.width; ++j) {
            auto c = GetImageColor(image, j, i);
            float sumc = (c.r + c.g + c.b) * c.a / 255.f;
            if (sumc != 0)
                field[i][j] = FloodCell{1.0f, Vector2Zero(), Vector3{(float)c.r, (float)c.g, (float)c.b}};
        }
    }
    return field;
}

void moveCell(const Vector2& imgSz, const FloodField& field, int fromX, int fromY, float coeff, FloodField& field2) {
    const Vector2& vel = field[fromY][fromX].vel;
    Vector2 to = Vector2{std::clamp(fromX + vel.x + 0.5f, 0 + 0.5f, imgSz.x - 0.5f), std::clamp(fromY + vel.y + 0.5f, 0 + 0.5f, imgSz.y - 0.5f)};
    int nj = int(floor(to.x)), ni = int(floor(to.y));
    Vector2 off = { (vel.x == 0 ? 0.0f : to.x - nj > 0.5f ? 0.5f : -0.5f), (vel.y == 0 ? 0.0f : to.y - ni > 0.5f ? 0.5f : -0.5f) };
    Vector2 a = { nj + 0.5f + off.x - abs(off.x), ni + 0.5f + off.y - abs(off.y)};
    Vector2 b = { nj + 0.5f + off.x - abs(off.x), ni + 0.5f + off.y + abs(off.y)};
    Vector2 c = { nj + 0.5f + off.x + abs(off.x), ni + 0.5f + off.y - abs(off.y)};
    Vector2 d = { nj + 0.5f + off.x + abs(off.x), ni + 0.5f + off.y + abs(off.y)};
    float aScore = 1.0f / std::max(Vector2Length(a - to), EPS);
    float bScore = 1.0f / std::max(Vector2Length(b - to), EPS);
    float cScore = 1.0f / std::max(Vector2Length(c - to), EPS);
    float dScore = 1.0f / std::max(Vector2Length(d - to), EPS);
    float sumScore = aScore + bScore + cScore + dScore;
    float aCoeff = coeff * aScore / sumScore;
    float bCoeff = coeff * bScore / sumScore;
    float cCoeff = coeff * cScore / sumScore;
    float dCoeff = coeff * dScore / sumScore;
    auto from = field[fromY][fromX];
    from.vel = to - Vector2{float(fromX) + 0.5f, float(fromY) + 0.5f};
    field2[std::clamp(int(floor(a.y)), 0, int(imgSz.y - 1))][std::clamp(int(floor(a.x)), 0, int(imgSz.x - 1))] += from * aCoeff;
    field2[std::clamp(int(floor(b.y)), 0, int(imgSz.y - 1))][std::clamp(int(floor(b.x)), 0, int(imgSz.x - 1))] += from * bCoeff;
    field2[std::clamp(int(floor(c.y)), 0, int(imgSz.y - 1))][std::clamp(int(floor(c.x)), 0, int(imgSz.x - 1))] += from * cCoeff;
    field2[std::clamp(int(floor(d.y)), 0, int(imgSz.y - 1))][std::clamp(int(floor(d.x)), 0, int(imgSz.x - 1))] += from * dCoeff;
}

void applyGravity(const Vector2& imgSz, FloodField& field) {
    for (int i = 0; i < imgSz.y; ++i) {
        for (int j = 0; j < imgSz.x; ++j) {
            FloodCell& fc = field[i][j];
            if (fc.mass > EPS) {
                fc.vel += GRAVITY;
            }
        }
    }
}

void applyPressure(const Vector2& imgSz, FloodField& field) {
    for (int i = 0; i < imgSz.y; ++i) {
        for (int j = 0; j < imgSz.x; ++j) {
            FloodCell& fc = field[i][j];
            if (fc.mass > EPS) {
                Vector2 v = Vector2Zero();
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        if (di == 0 && dj == 0)
                            continue;
                        int ii = i + di, jj = j + dj;
                        float m = 1.0f;
                        if (jj >= 0 && jj < imgSz.x && ii >= 0 && ii < imgSz.y) 
                            m = field[ii][jj].mass;
                        v += (1.0f - m) * Vector2Normalize(Vector2{float(dj), float(di)});
                    }
                }
                fc.vel += PRESSURE * v / 8.0f;
            }
        }
    }
}

void moveField(const Vector2& imgSz, FloodField& field) {
    FloodField field2(imgSz.y, std::vector<FloodCell>(imgSz.x, FloodCell{0.0f, Vector2Zero(), Vector3Zero()}));
    for (int i = 0; i < imgSz.y; ++i) {
        for (int j = 0; j < imgSz.x; ++j) {
            FloodCell& fc = field[i][j];
            if (fc.mass > EPS) {
                moveCell(imgSz, field, j, i, 1.0f, field2);
            }
        }
    }
    field = field2;
}

void drawField(Image& img, FloodField& field) {
    for (int i = 0; i < img.height; ++i) {
        for (int j = 0; j < img.width; ++j) {
            auto cv = field[i][j].color;
            unsigned char brg = (unsigned char)(255.0f * std::clamp(field[i][j].mass, 0.0f, 1.f));
            Color c = {(unsigned char)(cv.x), (unsigned char)(cv.y), (unsigned char)(cv.z), brg};
            //Color c = {brg, brg, brg, 255};
            //Color c = {(unsigned char)(2550.f * field[i][j].vel.x), (unsigned char)(2550.f * field[i][j].vel.y), brg, 255};
            ImageDrawPixel(&img, j, i, c);
        }
    }
}

int main() {

    FloodField field = makeField(INPUT_IMG);
    Vector2 imgSz = {(float)field[0].size(), (float)field.size()};

    InitWindow(imgSz.x * SCALE, imgSz.y * SCALE, "FLOOD");
    SetTargetFPS(60);

    Image img = GenImageColor(imgSz.x, imgSz.y, BLANK);
    Image imgCopy = ImageCopy(img);
    Texture2D tex = LoadTextureFromImage(img);

    Vector2 offset = Vector2Zero(), curOffset = Vector2Zero();
    float scale = 1.0f;
    Vector2 lastMouseGrabPos = GetMousePosition();

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R))
            field = makeField(INPUT_IMG);
        //if (IsKeyDown(KEY_SPACE)) {
            applyGravity(imgSz, field);
            applyPressure(imgSz, field);
            moveField(imgSz, field);
        //}
        drawField(img, field);

        UnloadImage(imgCopy);
        ImageCopy(img);
        Color *pixels = LoadImageColors(imgCopy);
        UpdateTexture(tex, pixels);         
        UnloadImageColors(pixels);              

        BeginDrawing();
        ClearBackground(BLACK);
        DrawTextureEx(tex, offset + curOffset, 0.0f, SCALE * scale, WHITE);
        EndDrawing();

        auto mpos = GetMousePosition();
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
            lastMouseGrabPos = mpos;
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
            offset = offset + curOffset;
        curOffset = IsMouseButtonDown(MOUSE_BUTTON_LEFT) ? (mpos - lastMouseGrabPos) : Vector2Zero();
        float newScale = scale + GetMouseWheelMove() * scale / 10.0f;
        offset -= ((mpos - offset) / scale) * (newScale - scale);
        scale = newScale;

        //Vector2 pix = (mpos - offset) / (scale * SCALE);
        //if (pix.x > 0 && pix.y > 0 && pix.x < imgSz.x && pix.y < imgSz.y) {
        //    auto& fc = field[int(floor(pix.y))][int(floor(pix.x))];
        //    std::cout << fc.mass << " v: " << fc.vel.x << " " << fc.vel.y << "\n";
        //}
    }

    CloseWindow();

    return 0;
}