#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"
#include "external/stb_image.h"

#include <cstddef>
#include <vector>

#define INPUT_IMG "res/input.png"

typedef struct FloodInput {
    Vector2 pos;
    float radius;
    unsigned int type;
} FloodInputCmd;

#define MAX_BUFFERED_INPUTS 64
typedef struct FloodInputSSBO {
    unsigned int count;
    FloodInput commands[MAX_BUFFERED_INPUTS];
} FloodInputSSBO;

struct FloodField {
    int W, H;
    unsigned int colSSBO[2], masVelSSBO[2], fixSSBO, transferSSBO;
    unsigned int transferProgram, logicProgram;
    FloodInputSSBO transferBuffer;
    Shader renderShader;
    Texture2D curTex;
    bool firstInitDone = false;
    int img_idx = 0;

    FloodField(const Image& img) {
        init(img);
    }

    ~FloodField() {
        for (int i = 0; i < 2; ++i) {
            rlUnloadShaderBuffer(colSSBO[i]);
            rlUnloadShaderBuffer(masVelSSBO[i]);
        }
        rlUnloadShaderBuffer(fixSSBO);
        rlUnloadShaderBuffer(transferSSBO);
        rlUnloadShaderProgram(transferProgram);
        rlUnloadShaderProgram(logicProgram);
        UnloadTexture(curTex);
        UnloadShader(renderShader);
    }

    void init(const Image& img, int w = 0, int h = 0) {
        if (!firstInitDone) {
            char *logicCode = LoadFileText("res/logic.glsl");
            unsigned int logicShader = rlCompileShader(logicCode, RL_COMPUTE_SHADER);
            logicProgram = rlLoadComputeShaderProgram(logicShader);

            renderShader = LoadShader(NULL, "res/render.glsl");

            char *transferCode = LoadFileText("res/transfer.glsl");
            unsigned int transferShader = rlCompileShader(transferCode, RL_COMPUTE_SHADER);
            transferProgram = rlLoadComputeShaderProgram(transferShader);
            UnloadFileText(transferCode);

            transferBuffer = { 0 };

            firstInitDone = true;
        } else {
            for (int i = 0; i < 2; ++i) {
                rlUnloadShaderBuffer(colSSBO[i]);
                rlUnloadShaderBuffer(masVelSSBO[i]);
             }
            rlUnloadShaderBuffer(fixSSBO);
        }
        bool fullscreen = (w != 0);
        auto sz = fullscreen ? Vector2{float(w), float(h)} : Vector2{float(img.width), float(img.height)};
        W = int(sz.x);
        H = int(sz.y);
        int LEFT = fullscreen ? int(sz.x * 0.5f - img.width * 0.5f) : 0;
        int TOP = fullscreen ? int(sz.y * 0.5f - img.height * 0.5f) : 0;
        SetShaderValue(renderShader, GetShaderLocation(renderShader, "RESOLUTION"), &sz, SHADER_UNIFORM_VEC2);

        std::vector<float> buffer(W * H * 3, 0.0f);
        for (int i = 0; i < img.height; ++i) {
            for (int j = 0; j < img.width; ++j) {
                auto c = GetImageColor(img, j, i);
                buffer[(TOP + i) * W * 3 + (LEFT + j) * 3 + 0] = c.r / 255.0f;
                buffer[(TOP + i) * W * 3 + (LEFT + j) * 3 + 1] = c.g / 255.0f;
                buffer[(TOP + i) * W * 3 + (LEFT + j) * 3 + 2] = c.b / 255.0f;
            }
        }
        colSSBO[0] = rlLoadShaderBuffer(W * H * 3 * sizeof(float), buffer.data(), RL_DYNAMIC_COPY);
        std::vector<unsigned char> buffer2(W * H, 0);
        for (int i = 0; i < img.height; ++i) {
            for (int j = 0; j < img.width; ++j) {
                auto c = GetImageColor(img, j, i);
                Vector2 vel = Vector2Zero();
                float mass = (c.r + c.g + c.b) / (3.0f * 255.0f);
                bool wall = c.r == 127 && c.g == 127 && c.b == 127;
                buffer[(TOP + i) * W * 3 + (LEFT + j) * 3 + 0] = vel.x;
                buffer[(TOP + i) * W * 3 + (LEFT + j) * 3 + 1] = vel.y;
                buffer[(TOP + i) * W * 3 + (LEFT + j) * 3 + 2] = mass;
                buffer2[(TOP + i) * W + (LEFT + j)] = wall;
            }
        }
        masVelSSBO[0] = rlLoadShaderBuffer(W * H * 3 * sizeof(float), buffer.data(), RL_DYNAMIC_COPY);
        colSSBO[1] = rlLoadShaderBuffer(W * H * 3 * sizeof(float), NULL, RL_DYNAMIC_COPY);
        masVelSSBO[1] = rlLoadShaderBuffer(W * H * 3 * sizeof(float), NULL, RL_DYNAMIC_COPY);
        fixSSBO = rlLoadShaderBuffer(W * H * sizeof(unsigned char), buffer2.data(), RL_DYNAMIC_COPY);

        auto tmp = GenImageColor(sz.x, sz.y, RED);
        curTex = LoadTextureFromImage(tmp);
        UnloadImage(tmp);
        img_idx = 0;
    }

    void swap() {
        img_idx = !img_idx;
    }

    void compute() {
        rlEnableShader(logicProgram);
        rlBindShaderBuffer(colSSBO[img_idx], 1);
        rlBindShaderBuffer(masVelSSBO[img_idx], 2);
        rlBindShaderBuffer(colSSBO[!img_idx], 3);
        rlBindShaderBuffer(masVelSSBO[!img_idx], 4);
        rlBindShaderBuffer(fixSSBO, 5);
        rlComputeShaderDispatch(W, H, 1);
        rlDisableShader();
    }

    void render() {
        rlBindShaderBuffer(colSSBO[img_idx], 1);
        BeginDrawing();
            ClearBackground(RED);
            BeginShaderMode(renderShader);
            DrawTexture(curTex, 0, 0, WHITE);
            EndShaderMode();
            //DrawRectangleLines(GetMouseX() - brushSize/2, GetMouseY() - brushSize/2, brushSize, brushSize, RED);
            //DrawText("Use Mouse wheel to increase/decrease brush size", 10, 10, 20, WHITE);
            //DrawFPS(GetScreenWidth() - 100, 10);
        EndDrawing();
    }
};

int main() {
    auto image = LoadImage(INPUT_IMG);
    InitWindow(image.width, image.height, "FLOOD1");
    SetTargetFPS(60);

    auto field = FloodField(image);

    bool gravityEnabled = true;
    bool mode = false;
    
    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R))
            field.init(image);

        //field.compute();
        //field.swap();
        field.render();

        if (IsKeyPressed(KEY_G))
            gravityEnabled = !gravityEnabled;
        if (IsKeyPressed(KEY_Q))
            mode = !mode;
        if (IsKeyPressed(KEY_F)) {
            if (!IsWindowFullscreen()) {
                field.init(image, GetMonitorWidth(0), GetMonitorHeight(0));
                SetWindowSize(field.W, field.H);
            }
            ToggleFullscreen();
            if (!IsWindowFullscreen()) {
                field.init(image);
                SetWindowSize(field.W, field.H);
            }
        }
    }

    UnloadImage(image);
    CloseWindow();

    return 0;
}