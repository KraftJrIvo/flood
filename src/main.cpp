#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "raylib.h"
#include "raymath.h"

#include "vec_ops.h"

#define INPUT_IMG "res/test1.png"
#define GRAVITY Vector2{0.0f, 0.01f}
#define EPS 1e-6
//#define RADIUS 1.0f

#define SCALE 4.0f

struct FloodCell {
    float mass;
    Vector2 vel;
    Vector3 color;
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

FloodCell interpolateCell(const Vector2& imgSz, const FloodField& field, const Vector2& pos) {
    FloodCell fc{0.0f, Vector2Zero(), Vector3Zero()};
    if (pos.x < 0 || pos.x >= imgSz.x || pos.y < 0 || pos.y >= imgSz.y)
        return fc;
    int xf = floor(pos.x), yf = floor(pos.y);
    int xc = ceil(pos.x), yc = ceil(pos.y);
    auto lt = field[yf][xf];
    auto rt = field[yf][xc];
    auto lb = field[yc][xf];
    auto rb = field[yc][xc];
    float tx = pos.x - xf, ty = pos.y - yf;
    fc.mass = (lt.mass + tx * (rt.mass - lt.mass)) + ty * ((lb.mass + tx * (rb.mass - lb.mass)) - (lt.mass + tx * (rt.mass - lt.mass)));
    fc.vel = (lt.vel + tx * (rt.vel - lt.vel)) + ty * ((lb.vel + tx * (rb.vel - lb.vel)) - (lt.vel + tx * (rt.vel - lt.vel)));
    fc.color = (lt.color + tx * (rt.color - lt.color)) + ty * ((lb.color + tx * (rb.color - lb.color)) - (lt.color + tx * (rt.color - lt.color)));
}

void updateField(const Vector2& imgSz, FloodField& field) {
    FloodField deltaField(imgSz.y, std::vector<FloodCell>(imgSz.x, FloodCell{0.0f, Vector2Zero(), Vector3Zero()}));
    for (int i = 0; i < imgSz.y; ++i) {
        for (int j = 0; j < imgSz.x; ++j) {
            FloodCell& fc = field[i][j];
            volatile FloodCell ffc = field[i][j];
            if (fc.mass > EPS) {
                //float masum = 0.0f;
                //for (int di = -1; di <= 1; ++di) {
                //    for (int dj = -1; dj <= 1; ++dj) {
                //        int ii = i + di, jj = j + dj;
                //        if (jj >= 0 && jj < imgSz.x && ii >= 0 && ii < imgSz.y)
                //            masum += field[i + di][j + dj].mass;
                //        else
                //            masum += 1.0f;
                //    }
                //}
                Vector2 v = Vector2Zero();
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        int ii = i + di, jj = j + dj;
                        if (jj >= 0 && jj < imgSz.x && ii >= 0 && ii < imgSz.y)
                            v += (1.0f - (field[i + di][j + dj].mass/* + fc.mass*/)) * Vector2{float(dj), float(di)};
                        else
                            v += -2.0f * Vector2{float(dj), float(di)};
                    }
                }
                fc.vel = fc.vel + v / 50.0f + GRAVITY;// - GRAVITY * (std::clamp(fc.mass - 1.0f, 0.0f, 2.0f));
                //fc.vel = Vector2{std::clamp(fc.vel.x, -0.75f, 0.75f), std::clamp(fc.vel.y, -0.75f, 0.75f)};
                Vector2 p = Vector2{j + 0.5f, i + 0.5f} + fc.vel;
                int nj = int(floor(p.x)), ni = int(floor(p.y));
                Vector2 off = { (p.x - nj > 0.5f ? 0.5f : -0.5f), (p.y - ni > 0.5f ? 0.5f : -0.5f) };
                Vector2 a = { nj + 0.5f + off.x - 0.5f, ni + 0.5f + off.y - 0.5f};
                Vector2 b = { nj + 0.5f + off.x - 0.5f, ni + 0.5f + off.y + 0.5f};
                Vector2 c = { nj + 0.5f + off.x + 0.5f, ni + 0.5f + off.y - 0.5f};
                Vector2 d = { nj + 0.5f + off.x + 0.5f, ni + 0.5f + off.y + 0.5f};
                float aDist = Vector2Length(a - p);
                float bDist = Vector2Length(b - p);
                float cDist = Vector2Length(c - p);
                float dDist = Vector2Length(d - p);
                float maxDist = std::max(std::max(aDist, bDist), std::max(cDist, dDist));
                float sumDist = 4.0f * maxDist - aDist - bDist - cDist - dDist;
                deltaField[i][j].mass -= fc.mass;
                deltaField[i][j].color -= fc.color;
                float dbg = 0.0f;
                int ax = int(floor(a.x)), ay = int(floor(a.y));
                if (ax >= 0 && ax < imgSz.x && ay >= 0 && ay < imgSz.y) {
                    deltaField[ay][ax].mass += fc.mass * (maxDist-aDist) / sumDist;
                    dbg += fc.mass * (maxDist-aDist) / sumDist;
                    deltaField[ay][ax].color = field[ay][ax].color + (fc.color - field[ay][ax].color) * (maxDist-aDist) / sumDist;
                } else {
                    deltaField[i][j].mass += fc.mass * (maxDist-aDist) / sumDist;
                    field[i][j].vel = field[i][j].vel * 0.1f;
                }
                int bx = int(floor(b.x)), by = int(floor(b.y));
                if (bx >= 0 && bx < imgSz.x && by >= 0 && by < imgSz.y) {
                    deltaField[by][bx].mass += fc.mass * (maxDist - bDist) / sumDist;
                    dbg += fc.mass * (maxDist - bDist) / sumDist;
                    deltaField[by][bx].color = field[by][bx].color + (fc.color - field[by][bx].color) * (maxDist - bDist) / sumDist;
                } else {
                    deltaField[i][j].mass += fc.mass * (maxDist - bDist) / sumDist;
                    field[i][j].vel = field[i][j].vel * 0.1f;
                }
                int cx = int(floor(c.x)), cy = int(floor(c.y));
                if (cx >= 0 && cx < imgSz.x && cy >= 0 && cy < imgSz.y) {
                    deltaField[cy][cx].mass += fc.mass * (maxDist - cDist) / sumDist;
                    dbg += fc.mass * (maxDist - cDist) / sumDist;
                    deltaField[cy][cx].color = field[cy][cx].color + (fc.color - field[cy][cx].color) * (maxDist - cDist) / sumDist;
                } else {
                    deltaField[i][j].mass += fc.mass * (maxDist - cDist) / sumDist;
                    field[i][j].vel = field[i][j].vel * 0.1f;
                }
                int dx = int(floor(d.x)), dy = int(floor(d.y));
                if (dx >= 0 && dx < imgSz.x && dy >= 0 && dy < imgSz.y) {
                    deltaField[dy][dx].mass += fc.mass * (maxDist - dDist) / sumDist;
                    dbg += fc.mass * (maxDist - dDist) / sumDist;
                    deltaField[dy][dx].color = field[dy][dx].color + (fc.color - field[dy][dx].color) * (maxDist - dDist) / sumDist;
                } else {
                    deltaField[i][j].mass += fc.mass * (maxDist - dDist) / sumDist;
                    field[i][j].vel = field[i][j].vel * 0.1f;
                }
                //std::cout << field[ay][ax].mass << " " << deltaField[ay][ax].mass << "\n";
                //std::cout << field[by][bx].mass << " " << deltaField[by][bx].mass << "\n";
                //std::cout << field[cy][cx].mass << " " << deltaField[cy][cx].mass << "\n";
                //std::cout << field[dy][dx].mass << " " << deltaField[dy][dx].mass << "\n";
                //std::cout << "";
            }
        }
    }
    for (int i = 0; i < imgSz.y; ++i) {
        for (int j = 0; j < imgSz.x; ++j) {
            field[i][j].mass += deltaField[i][j].mass;
            field[i][j].color += deltaField[i][j].color;
        }
    }
}

void drawField(Image& img, FloodField& field) {
    for (int i = 0; i < img.height; ++i) {
        for (int j = 0; j < img.width; ++j) {
            auto cv = field[i][j].color;
            //Color c = {(unsigned char)(cv.x), (unsigned char)(cv.y), (unsigned char)(cv.z), 255};
            unsigned char brg = (unsigned char)(255.0f * std::clamp(field[i][j].mass, 0.0f, 1.f));
            Color c = {brg, brg, brg, 255};
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
        //if (IsKeyPressed(KEY_SPACE))
            updateField(imgSz, field);
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

        Vector2 pix = (mpos - offset) / (scale * SCALE);
        if (pix.x > 0 && pix.y > 0 && pix.x < imgSz.x && pix.y < imgSz.y) {
            auto& fc = field[int(floor(pix.y))][int(floor(pix.x))];
            std::cout << fc.mass << " " << fc.color.x << " " << fc.color.y << " " << fc.color.z << "\n";
        }
    }

    CloseWindow();

    return 0;
}