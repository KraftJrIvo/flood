#version 430

layout (local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

layout(std430, binding = 1) readonly restrict buffer colorsLaySrc {
    uint colorsBufSrc[];
};

layout(std430, binding = 2) readonly restrict buffer masVelLaySrc {
    uint masVelBufSrc[];
};

layout(std430, binding = 3) writeonly restrict buffer colorsLayDst {
    uint colorsDst[];
};

layout(std430, binding = 4) writeonly restrict buffer masVelLayDSt {
    uint masVelDst[];
};

layout(std430, binding = 5) readonly restrict buffer fixLaySrc {
    uint fixBufSrc[];
};

void main()
{
    //uint neighbourCount = 0;
    //uint x = gl_GlobalInvocationID.x;
    //uint y = gl_GlobalInvocationID.y;
//
    //neighbourCount += fetchGol(x - 1, y - 1);   // Top left
    //neighbourCount += fetchGol(x, y - 1);       // Top middle
    //neighbourCount += fetchGol(x + 1, y - 1);   // Top right
    //neighbourCount += fetchGol(x - 1, y);       // Left
    //neighbourCount += fetchGol(x + 1, y);       // Right
    //neighbourCount += fetchGol(x - 1, y + 1);   // Bottom left
    //neighbourCount += fetchGol(x, y + 1);       // Bottom middle   
    //neighbourCount += fetchGol(x + 1, y + 1);   // Bottom right
//
    //if (neighbourCount == 3) setGol(x, y, 1);
    //else if (neighbourCount == 2) setGol(x, y, fetchGol(x, y));
    //else setGol(x, y, 0);
}