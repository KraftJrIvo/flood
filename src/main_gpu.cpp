#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"
#include "external/stb_image.h"

#include "vec_ops.h"
#include <string>
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

//Image loadImageRGBAfloat32(const Image& imgsrc, int w, int h, const Vector2& topleft) {        
//        int channels = 4;
//        int len = channels * w * h;
//        std::vector<float> buffer(len, 0.0f);
//        int left = topleft.x;
//        int top = topleft.y;
//        for (int i = 0; i < imgsrc.height; ++i) {
//            for (int j = 0; j < imgsrc.width; ++j) {
//                auto c = GetImageColor(imgsrc, j, i);
//                Vector2 vel = Vector2Zero();
//                float mass = (c.r + c.g + c.b) / (3.0f * 255.0f);
//                bool wall = c.r == 127 && c.g == 127 && c.b == 127;
//                buffer[(top + i) * w * channels + (left + j) * channels + 0] = vel.x;
//                buffer[(top + i) * w * channels + (left + j) * channels + 1] = vel.y;
//                buffer[(top + i) * w * channels + (left + j) * channels + 2] = mass;
//                buffer[(top + i) * w * channels + (left + j) * channels + 3] = wall;
//            }
//        }
//        Image img;
//        img.data = stbi_loadf_from_memory((unsigned char*)buffer.data(), len * sizeof(float), &w, &h, &channels, channels);
//        img.mipmaps = 1;
//        img.format = PIXELFORMAT_UNCOMPRESSED_R32G32B32A32;
//        img.width = w;
//        img.height = h;
//        return img;
//}
//
//struct FloodField {
//    Vector2 sz;
//    //Shader shader;
//    //Texture2D colors[2];
//    //Texture2D massVelFix[2];
//    //Texture2D curTex;
//    //int img_idx = 0;
//    bool shaderLoaded = false;
//
//    unsigned int colSSBO, masVelSSBO, fixSSBO;
//
//    FloodField(const Image& img) {
//        init(img);
//    }
//
//    void init(const Image& img, const Vector2& winSz = Vector2Zero()) {
//        if (!shaderLoaded) {
//            //shader = LoadShader(nullptr, "res/flood.frag");
//            shaderLoaded = true;
//        } else {
//            rlUnloadShaderBuffer(ssboA);
//            rlUnloadShaderBuffer(ssboB);
//            rlUnloadShaderBuffer(ssboTransfert);
//
//            // Unload compute shader programs
//            rlUnloadShaderProgram(golTransfertProgram);
//            rlUnloadShaderProgram(golLogicProgram);
//
//            UnloadTexture(whiteTex);            // Unload white texture
//            UnloadShader(golRenderShader);      // Unload rendering fragment shader 
//        }
//        bool fullscreen = (winSz.x != 0);
//        sz = fullscreen ? winSz : Vector2{float(img.width), float(img.height)};
//        Image colimg = (fullscreen) ? GenImageColor(sz.x, sz.y, BLACK) : img;
//        Vector2 topleft = Vector2Zero();
//        if (fullscreen) {
//            topleft = Vector2{float(int(winSz.x * 0.5f - img.width * 0.5f)), float(int(winSz.y * 0.5f - img.height * 0.5f))};
//            ImageDraw(&colimg, img, Rectangle{0, 0, float(img.width), float(img.height)}, Rectangle{ topleft.x, topleft.y, float(img.width), float(img.height) }, WHITE);
//        }
//        Image massVelFixImg = loadImageRGBAfloat32(img, sz.x, sz.y, topleft);        
//        for (int i = 0; i < 2; ++i) {
//            colors[i] = LoadTextureFromImage(colimg);
//            massVelFix[i] = LoadTextureFromImage(massVelFixImg);
//            //SetShaderValueTexture(shader, GetShaderLocation(shader, ("colors_" + std::to_string(i)).c_str()), colors[i]);
//            //SetShaderValueTexture(shader, GetShaderLocation(shader, ("masVelFix_" + std::to_string(i)).c_str()), massVelFix[i]);
//            rlBindImageTexture(3 + i * 2 + 0, colors[i].id, RL_PIXELFORMAT_UNCOMPRESSED_R8G8B8A8, false);
//            rlBindImageTexture(3 + i * 2 + 1, massVelFix[i].id, RL_PIXELFORMAT_UNCOMPRESSED_R32G32B32A32, false);
//        }
//        curTex = LoadTextureFromImage(GenImageColor(sz.x, sz.y, BLACK));
//        SetShaderValue(shader, GetShaderLocation(shader, "RESOLUTION"), &sz, SHADER_ATTRIB_VEC2);
//        img_idx = 0;
//        swap(true);
//    }
//
//    Texture2D getColorTex() {
//        return colors[img_idx];
//    }
//
//    Texture2D getMassVelFixTex() {
//        return massVelFix[img_idx];
//    }
//
//    void swap(bool first = false) {
//        if (!first) img_idx = !img_idx;
//        SetShaderValue(shader, GetShaderLocation(shader, "IMG_IDX"), &img_idx, SHADER_UNIFORM_INT);
//    }
//};

struct FloodField {
    Vector2 sz;
    unsigned int colSSBO, masVelSSBO, fixSSBO, transferSSBO;
    unsigned int transferProgram, logicProgram;
    FloodInputSSBO transferBuffer;
    Shader renderShader;
    Texture2D curTex;
    bool firstInitDone = false;

    FloodField(const Image& img) {
        init(img);
    }

    ~FloodField() {
        rlUnloadShaderBuffer(colSSBO);
        rlUnloadShaderBuffer(masVelSSBO);
        rlUnloadShaderBuffer(fixSSBO);
        rlUnloadShaderBuffer(transferSSBO);
        rlUnloadShaderProgram(transferProgram);
        rlUnloadShaderProgram(logicProgram);
        UnloadTexture(curTex);
        UnloadShader(renderShader);
    }

    void init(const Image& img, const Vector2& winSz = Vector2Zero()) {
        if (!firstInitDone) {
            char *logicCode = LoadFileText("res/logic.glsl");
            unsigned int logicShader = rlCompileShader(logicCode, RL_COMPUTE_SHADER);
            logicProgram = rlLoadComputeShaderProgram(logicShader);

            renderShader = LoadShader(NULL, "res/render.glsl");

            char *transferCode = LoadFileText("resources/shaders/glsl430/gol_transfert.glsl");
            unsigned int transferShader = rlCompileShader(transferCode, RL_COMPUTE_SHADER);
            transferProgram = rlLoadComputeShaderProgram(transferShader);
            UnloadFileText(transferCode);

            transferBuffer = { 0 };

            firstInitDone = true;
        } else {
            rlUnloadShaderBuffer(colSSBO);
            rlUnloadShaderBuffer(masVelSSBO);
            rlUnloadShaderBuffer(fixSSBO);
        }
        bool fullscreen = (winSz.x != 0);
        sz = fullscreen ? winSz : Vector2{float(img.width), float(img.height)};
        int W = int(sz.x), H = int(sz.y);
        colSSBO = rlLoadShaderBuffer(W * H * 3 * sizeof(float), NULL, RL_DYNAMIC_COPY);
        masVelSSBO = rlLoadShaderBuffer(W * H * 3 * sizeof(float), NULL, RL_DYNAMIC_COPY);
        fixSSBO = rlLoadShaderBuffer(W * H * sizeof(unsigned char), NULL, RL_DYNAMIC_COPY);
    }
};

int main() {
    auto image = LoadImage(INPUT_IMG);
    InitWindow(image.width, image.height, "FLOOD");
    SetTargetFPS(60);

    auto field = FloodField(image);

    bool gravityEnabled = true;
    bool mode = false;
    
    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_R))
            field.init(image);

        BeginDrawing();
        ClearBackground(BLACK);
        BeginShaderMode(field.shader);
        
        //SetShaderValueTexture(field.shader, GetShaderLocation(field.shader, "colors_0"), field.colors[0]);
        //SetShaderValueTexture(field.shader, GetShaderLocation(field.shader, "masVelFix_0"), field.massVelFix[0]);
        //rlBindImageTexture(;// (field.colors[0].id, GetShaderLocation(field.shader, "colors_0"), false)
        DrawTexture(field.curTex, 0, 0, WHITE);
        EndShaderMode();
        EndDrawing();

        field.swap();

        if (IsKeyPressed(KEY_G))
            gravityEnabled = !gravityEnabled;
        if (IsKeyPressed(KEY_Q))
            mode = !mode;
        if (IsKeyPressed(KEY_F)) {
            if (!IsWindowFullscreen()) {
                field.init(image, {(float)GetMonitorWidth(0), (float)GetMonitorHeight(0)});
                SetWindowSize(field.sz.x, field.sz.y);
            }
            ToggleFullscreen();
            if (!IsWindowFullscreen()) {
                field.init(image);
                SetWindowSize(field.sz.x, field.sz.y);
            }
        }
    }

    CloseWindow();

    return 0;
}