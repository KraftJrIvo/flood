#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <array>
#include <algorithm>

#include "raylib.h"
#include "raymath.h"

#include "vec_ops.h"

#define INPUT_IMG "res/input.png"

#define EPS 1e-3f

#define GRAVITY Vector2{0.0f, 0.1f}

#define VEL_FADE 0.995f
#define PRESSURE 1.0f
#define MIN_VISC 0.1f
#define MAX_VISC 0.85f

#define M_RADIUS 30.0f
#define M_FORCE 0.3f

struct Matter {
    std::string name;
    uint8_t idx;
    Vector3 dftColor;
    float density;
    float viscosity;
    float compressability;
    float granularity;
};
typedef std::shared_ptr<Matter> MatterPtr;

struct MatterPortion {
    MatterPtr matter = nullptr;
    float mass = 0.0f;
    Vector2 vel = Vector2Zero();
    Vector3 color = Vector3Zero();
    void operator+=(const MatterPortion& b) {
        matter = b.matter;
        vel = VEL_FADE * (vel * mass / std::max(mass + b.mass, EPS) + b.vel * b.mass / std::max(mass + b.mass, EPS));
        color = color * mass / std::max(mass + b.mass, EPS) + b.color * b.mass / std::max(mass + b.mass, EPS);
        mass += b.mass;
    }
    MatterPortion operator*(float coeff) const {
        MatterPortion mc;
        mc.matter = matter;
        mc.vel = vel;
        mc.mass = mass * coeff;
        mc.color = color;
        return mc;
    }
};

template <int N>
struct MatterCell {
    std::array<MatterPortion, N> content;
    bool fixed = false;
    float mass() const {
        float mass = 0.0f;
        for (auto& mp : content)
            mass += mp.mass;
        return mass;
    }
    Vector2 vel() const {
        Vector2 vel = Vector2Zero();
        for (auto& mp : content)
            vel += mp.vel;
        return vel;
    }
};

template <int N>
struct MatterField {
    Vector2 sz;
    std::vector<std::vector<MatterCell<N>>> cells;
};

template <int N>
MatterField<N> makeField(std::string imagePath, const std::vector<MatterPtr>& matters) {
    auto image = LoadImage(imagePath.c_str());
    MatterField<N> field{Vector2{float(image.width), float(image.height)}, {(size_t)image.height, std::vector<MatterCell<N>>{(size_t)image.width, MatterCell<N>()}}};
    for (int i = 0; i < image.height; ++i) {
        for (int j = 0; j < image.width; ++j) {
            auto c = GetImageColor(image, j, i);
            float sumc = (c.r + c.g + c.b) * c.a / 255.f;
            bool gray = ((c.r == 127) && (c.g == 127) && (c.b == 127));
            if (sumc != 0) {
                auto mid = gray ? 1 : 0;
                field.cells[i][j].content[mid] = MatterPortion{matters[mid], matters[mid]->density, Vector2Zero(), Vector3{(float)c.r, (float)c.g, (float)c.b}};
                field.cells[i][j].fixed = gray;
            }
        }
    }
    return field;
}

std::pair<MatterPortion, MatterPortion> splitMatterPortion(const MatterPortion& mp, float ratio) {
    auto m = int(mp.mass / mp.matter->granularity);
    auto n = int(ratio * mp.mass / mp.matter->granularity);
    MatterPortion a = mp; a.mass = n * mp.matter->granularity;
    MatterPortion b = mp; b.mass = (m - n) * mp.matter->granularity;
    return {a, b};
}

template <int N>
void walkAndPlace(const MatterField<N>& field, MatterField<N>& field2, MatterPortion mp, const Vector2& pos, const Vector2& vel, const Vector2& dir, int fx, int fy) {
    Vector2 end = pos;
    Vector2 step = {(dir.x > 0) ? 1.0f : -1.0f, (dir.y > 0) ? 1.0f : -1.0f};
    int j = int(floor(pos.x)), i = int(floor(pos.y));
    while (true) {
        if (j == fx && i == fy) {
            mp.vel = vel;
            field2.cells[i][j].content[mp.matter->idx] += mp;
            break;
        }
        Vector2 targ = Vector2{(dir.x > 0) ? ceil(end.x) : floor(end.x), (dir.y > 0) ? ceil(end.y) : floor(end.y)} - end;
        Vector2 targX = {targ.x, targ.x * dir.y/(dir.x+EPS)};
        Vector2 targY = {targ.y * dir.x/(dir.y+EPS), targ.y};
        Vector2 tmax = {Vector2Length(targX), Vector2Length(targY)};
        Vector2 newend = end + Vector2{(tmax.x < tmax.y) ? targX.x : targY.x, (tmax.x < tmax.y) ? targX.y : targY.y} + step * Vector2{EPS, EPS};
        int nj = int(floor(newend.x)), ni = int(floor(newend.y));
        if (field.cells[ni][nj].fixed/* || (field2.cells[ni][nj].content[mp.matter->idx].mass * mp.matter->compressability) > mp.matter->density*/) {
            mp.vel = vel - Vector2{float(fx - j), float(fy - i)};
            field2.cells[i][j].content[mp.matter->idx] += mp;
            break;
        }
        i = ni; j = nj;
        end = newend;
    }
}

template <int N>
void moveCell(const Vector2& imgSz, const MatterField<N>& field, int fromX, int fromY, MatterField<N>& field2) {
    for (int k = 0; k < N; ++k) {
        auto from = field.cells[fromY][fromX].content[k];    
        if (from.matter) {
            const Vector2& vel = from.vel;
            Vector2 pos = Vector2{float(fromX) + 0.5f, float(fromY) + 0.5f};
            Vector2 to = Vector2{std::clamp(fromX + vel.x + 0.5f, 0 + 0.5f, imgSz.x - 0.5f), std::clamp(fromY + vel.y + 0.5f, 0 + 0.5f, imgSz.y - 0.5f)};
            int nj = int(floor(to.x)), ni = int(floor(to.y));
            Vector2 off = { (vel.x == 0 ? 0.0f : to.x - nj > 0.5f ? 0.5f : -0.5f), (vel.y == 0 ? 0.0f : to.y - ni > 0.5f ? 0.5f : -0.5f) };
            Vector2 a = { std::clamp(nj + 0.5f + off.x - abs(off.x), 0 + 0.5f, imgSz.x - 0.5f), std::clamp(ni + 0.5f + off.y - abs(off.y), 0 + 0.5f, imgSz.y - 0.5f)};
            Vector2 b = { std::clamp(nj + 0.5f + off.x - abs(off.x), 0 + 0.5f, imgSz.x - 0.5f), std::clamp(ni + 0.5f + off.y + abs(off.y), 0 + 0.5f, imgSz.y - 0.5f)};
            Vector2 c = { std::clamp(nj + 0.5f + off.x + abs(off.x), 0 + 0.5f, imgSz.x - 0.5f), std::clamp(ni + 0.5f + off.y - abs(off.y), 0 + 0.5f, imgSz.y - 0.5f)};
            Vector2 d = { std::clamp(nj + 0.5f + off.x + abs(off.x), 0 + 0.5f, imgSz.x - 0.5f), std::clamp(ni + 0.5f + off.y + abs(off.y), 0 + 0.5f, imgSz.y - 0.5f)};
            std::vector<Vector2> pts = {a, b, c, d}; 
            float score[5] = {0,0,0,0,0};
            for (int i = 0; i < 4; ++i) {
                score[i] = 1.0f / std::max(Vector2Length(pts[i] - to), EPS);
                score[4] += score[i];
            }
            for (int i = 0; i < 4; ++i) {
                ni = std::clamp(int(floor(pts[i].y)), 0, int(imgSz.y - 1));
                nj = std::clamp(int(floor(pts[i].x)), 0, int(imgSz.x - 1));
                float ratio = score[i] / score[4];
                if (from.matter->density - field2.cells[ni][nj].mass() > from.mass) {
                    from.vel = to - pos;
                    field2.cells[ni][nj].content[from.matter->idx] += from * ratio;
                } else {
                    walkAndPlace(field, field2, from * ratio/*splitMatterPortion(from, ratio).first*/, pos, to - pos, Vector2Normalize(pts[i] - pos), nj, ni);
                }
            }
        }
    }
}

template <int N>
void applyGravity(MatterField<N>& field) {
    for (int i = 0; i < field.sz.y; ++i) {
        for (int j = 0; j < field.sz.x; ++j) {
            if (field.cells[i][j].fixed)
                continue;
            for (int k = 0; k < N; ++k) {
                auto& mp = field.cells[i][j].content[k];
                if (mp.mass > EPS)
                    mp.vel += GRAVITY;
            }
        }
    }
}

template <int N>
void applyMouse(MatterField<N>& field, Vector2 pos, float radius, float force) {
    for (int i = 0; i < field.sz.y; ++i) {
        for (int j = 0; j < field.sz.x; ++j) {
            Vector2 ji = Vector2{j + 0.5f, i + 0.5f};
            float dist = Vector2Distance(ji, pos);
            if (field.cells[i][j].fixed || dist > radius)
                continue;
            for (int k = 0; k < N; ++k) {
                auto& mp = field.cells[i][j].content[k];
                if (mp.mass > EPS)
                    mp.vel += (1.0f - (dist / radius)) * Vector2Normalize(pos - ji) * force;
            }
        }
    }
}

template <int N>
void applyPressure(MatterField<N>& field) {
    for (int i = 0; i < field.sz.y; ++i) {
        for (int j = 0; j < field.sz.x; ++j) {
            if (field.cells[i][j].fixed)
                continue;
            for (int k = 0; k < N; ++k) {
                MatterPortion& mp = field.cells[i][j].content[k];
                if (mp.matter) {
                    float maxMass = 0.0f;
                    float nMass = 0.0f;
                    Vector2 v = Vector2Zero();
                    Vector2 v2 = Vector2Zero();
                    for (int di = -1; di <= 1; ++di) {
                        for (int dj = -1; dj <= 1; ++dj) {
                            int ii = i + di, jj = j + dj;
                            float mass2 = mp.matter->density * 3.0f;
                            if (jj >= 0 && jj < field.sz.x && ii >= 0 && ii < field.sz.y) {
                                if (field.cells[ii][jj].fixed) {
                                    mass2 = mp.matter->density;
                                } else {
                                    mass2 = std::clamp(field.cells[ii][jj].mass(), 0.0f, mp.matter->density * 50.0f);
                                    float massMP = field.cells[ii][jj].content[k].mass;
                                    if (massMP > EPS) {
                                        maxMass = std::max(maxMass, mass2);
                                        nMass += 1.0f;
                                        v2 += field.cells[ii][jj].vel();
                                    }
                                }
                            }
                            if (di == 0 && dj == 0)
                                continue;
                            if (mp.mass < mp.matter->density || mp.mass > mp.matter->density)
                                v += -sqrt(fabs(mp.matter->density - (mass2 + (mp.mass - mp.matter->density)))) * Vector2Normalize(Vector2{float(dj), float(di)});
                            else
                                v += -fabs(mp.matter->density - (mass2 + (mp.mass - mp.matter->density))) * Vector2Normalize(Vector2{float(dj), float(di)});
                        }
                    }
                    Vector2 vdelta = Vector2Zero();
                    vdelta += (PRESSURE + (MAX_VISC - mp.matter->viscosity)) * v / 8.0f;
                    vdelta += mp.matter->viscosity * (v2 * (1.0f / std::max(nMass, 1.0f)) - mp.vel);
                    mp.vel += vdelta;
                    for (int di = -1; di <= 1; ++di) {
                        for (int dj = -1; dj <= 1; ++dj) {
                            int ii = i + di, jj = j + dj;
                            if (jj >= 0 && jj < field.sz.x && ii >= 0 && ii < field.sz.y && !field.cells[ii][jj].fixed) {
                                field.cells[ii][jj].content[k].vel -= mp.matter->viscosity * (v2 * (1.0f / std::max(nMass, 1.0f)) - mp.vel) * 0.5f;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int N>
void moveField(MatterField<N>& field) {
    MatterField<N> field2{Vector2{float(field.sz.x), float(field.sz.y)}, {(size_t)field.sz.y, std::vector<MatterCell<N>>{(size_t)field.sz.x, MatterCell<N>()}}};
    for (int i = 0; i < field.sz.y; ++i) {
        for (int j = 0; j < field.sz.x; ++j) {
            MatterCell<N>& mc = field.cells[i][j];
            if (mc.fixed) {
                field2.cells[i][j] = mc;
            }
        }
    }
    for (int i = 0; i < field.sz.y; ++i) {
        for (int j = 0; j < field.sz.x; ++j) {
            MatterCell<N>& mc = field.cells[i][j];
            if (!mc.fixed && mc.mass() > EPS) {
                moveCell(field.sz, field, j, i, field2);
            }
        }
    }
    field = field2;
}

template <int N>
void drawField(Image& img, MatterField<N>& field) {
    ImageClearBackground(&img, BLACK);
    for (int i = 0; i < img.height; ++i) {
        for (int j = 0; j < img.width; ++j) {
            Vector3 cv = Vector3Zero();
            unsigned char brg = 0;
            float mass = 0;
            float vel = 0;
            for (int k = 0; k < N; ++k) {
                auto& mp = field.cells[i][j].content[k];
                if (mp.matter) {
                    brg = (unsigned char)std::clamp((brg + 255.0f * std::clamp(mp.mass, 0.0f, mp.matter->density) / mp.matter->density), float(brg), 255.0f);                    
                    Vector3 clr = mp.color;
                    cv = cv * mass / std::max(mass + mp.mass, EPS) + clr * mp.mass / std::max(mass + mp.mass, EPS);
                    mass += std::max((mp.mass - mp.matter->density) * PRESSURE, 0.0f);
                    vel += Vector2Length(mp.vel);
                }
            }
            //Color c = {(unsigned char)(cv.x), (unsigned char)(cv.y), (unsigned char)(cv.z), brg};
            Color c = {(unsigned char)(std::clamp(cv.x + vel * 30.0f + mass * 15.0f, 0.0f, 255.0f)), (unsigned char)(std::clamp(cv.y + vel * 30.0f, 0.0f, 255.0f)), (unsigned char)(std::clamp(cv.z + vel * 30.0f, 0.0f, 255.0f)), brg};
            ImageDrawPixel(&img, j, i, c);
        }
    }
}

int main() {

    std::vector<MatterPtr> matters = {
        std::make_shared<Matter>(Matter{"LIQUID", 0, Vector3{0.0f, 0.0f, 1.0f}, 1.0f, 1.0f, 0.5, EPS}),
        std::make_shared<Matter>(Matter{"WALL",   1, Vector3{0.5f, 0.5f, 0.5f}, 1.5f, 0.1f, 0, EPS})
    };
    auto field = makeField<2>(INPUT_IMG, matters);

    InitWindow(1, 1, "FLOOD");
    float SCALE = floor(GetMonitorHeight(0) / field.sz.y);

    SetWindowSize(field.sz.x * SCALE, field.sz.y * SCALE);
    SetWindowPosition((GetMonitorWidth(0) - field.sz.x * SCALE) * 0.5, (GetMonitorHeight(0) - field.sz.y * SCALE) * 0.5);
    SetTargetFPS(60);

    Image img = GenImageColor(field.sz.x, field.sz.y, BLANK);
    Texture2D tex = LoadTextureFromImage(img);

    Vector2 offset = Vector2Zero(), curOffset = Vector2Zero();
    float scale = 1.0f;
    //Vector2 lastMouseGrabPos = GetMousePosition();

    bool gravityEnabled = true;
    bool mode = false;
    
    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R))
            field = makeField<2>(INPUT_IMG, matters);
        if (gravityEnabled) applyGravity(field);
        applyMouse(field, (GetMousePosition() - offset) / (scale * SCALE), M_RADIUS, IsMouseButtonDown(1) ? M_FORCE : IsMouseButtonDown(0) ? -M_FORCE : 0.0f );
        applyPressure(field);
        moveField(field);
        drawField(img, field);

        Color *pixels = LoadImageColors(img);
        UpdateTexture(tex, pixels);         
        UnloadImageColors(pixels);     

        BeginDrawing();
        ClearBackground(BLACK);
        DrawTextureEx(tex, offset + curOffset, 0.0f, SCALE * scale, WHITE);
        EndDrawing();

        if (IsKeyPressed(KEY_G))
            gravityEnabled = !gravityEnabled;
        if (IsKeyPressed(KEY_Q))
            mode = !mode;
        matters[0]->viscosity = (mode ? MIN_VISC : MAX_VISC);

        // PAN + ZOOM
        //auto mpos = GetMousePosition();
        //if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
        //    lastMouseGrabPos = mpos;
        //if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
        //    offset = offset + curOffset;
        //curOffset = IsMouseButtonDown(MOUSE_BUTTON_LEFT) ? (mpos - lastMouseGrabPos) : Vector2Zero();
        //float newScale = scale + GetMouseWheelMove() * scale / 10.0f;
        //offset -= ((mpos - offset) / scale) * (newScale - scale);
        //scale = newScale;

        // DBG
        //Vector2 pix = (mpos - offset) / (scale * SCALE);
        //if (pix.x > 0 && pix.y > 0 && pix.x < field.sz.x && pix.y < field.sz.y) {
        //    auto& fc = field[int(floor(pix.y))][int(floor(pix.x))];
        //    std::cout << fc.mass << " v: " << fc.vel.x << " " << fc.vel.y << "\n";
        //}
    }

    CloseWindow();

    return 0;
}