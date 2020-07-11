#include <student/gpu.hpp>
#include <algorithm>


GPU::GPU()
{
}

BufferID GPU::createBuffer(uint64_t size)
{ 
    Buffer buffer{ size };
    BufferID id = buffer.GetID();
    buffers.emplace(id, std::move(buffer));
    return id; 
}
BufferID GPU::createBuffer(Buffer** _out, uint64_t size)
{
    Buffer buffer(size);
    BufferID id = buffer.GetID();
    buffers.emplace(id, std::move(buffer));
    *_out = &buffers.at(id); // std::c++11 pointers to map are not invalidated
    return id;
}
void GPU::deleteBuffer(BufferID buffer) 
{
    if (auto bufferit = buffers.find(buffer); bufferit != buffers.end())
    {
        buffers.erase(bufferit);
    }
}
void GPU::setBufferData(BufferID buffer, uint64_t offset, uint64_t size, void const* data) 
{
    if (auto bufferit = buffers.find(buffer); bufferit != buffers.end())
    {
        std::copy((uint8_t*)data, (uint8_t*)data + size, bufferit->second.GetRawData() + offset);
    }
}
void GPU::getBufferData(BufferID buffer, uint64_t offset, uint64_t size, void* data)
{
    if (auto bufferit = buffers.find(buffer); bufferit != buffers.end() && data)
    {
        auto pfirst = bufferit->second.GetRawData() + offset;
        std::copy(pfirst, pfirst + size, (uint8_t*)data);
    }
}
bool GPU::isBuffer(BufferID buffer) 
{
    return buffers.find(buffer) != buffers.end();
}


/// @}

ObjectID GPU::createVertexPuller()
{
    auto ptr = std::make_unique<InputLayout>();
    VertexPullerID ID = reinterpret_cast<VertexPullerID>(ptr.get());
    InputAssemblers.emplace(ID, std::move(ptr));
    return ID;
}
ObjectID GPU::createVertexPuller(InputLayout** _out)
{
    auto ptr = std::make_unique<InputLayout>();
    *_out = ptr.get();
    VertexPullerID ID = reinterpret_cast<VertexPullerID>(ptr.get());
    InputAssemblers.emplace(ID, std::move(ptr));
    return ID;
}
void GPU::deleteVertexPuller(VertexPullerID vao) 
{
    if (auto bufferit = InputAssemblers.find(vao); bufferit != InputAssemblers.end())
    {
        InputAssemblers.erase(bufferit);
    }
}
void GPU::setVertexPullerHead(VertexPullerID vao, uint32_t head, AttributeType type, uint64_t stride, uint64_t offset, BufferID buffer)
{
    if (auto bufferit = InputAssemblers.find(vao); bufferit != InputAssemblers.end())
    {
        bufferit->second->SetDescriptor(head, type, stride, offset, &buffers.at(buffer));//here
    }
}
void GPU::setVertexPullerIndexing(VertexPullerID vao, IndexType type, BufferID buffer) 
{
    if (auto bufferit = InputAssemblers.find(vao); bufferit != InputAssemblers.end())
    {
        bufferit->second->SetIndexing(&buffers.at(buffer), type);
    }
}
void GPU::enableVertexPullerHead(VertexPullerID vao, uint32_t head)
{
    if (auto bufferit = InputAssemblers.find(vao); bufferit != InputAssemblers.end())
    {
        bufferit->second->EnableHead(head);
    }
}
void GPU::disableVertexPullerHead(VertexPullerID vao, uint32_t head)
{
    if (auto bufferit = InputAssemblers.find(vao); bufferit != InputAssemblers.end())
    {
        bufferit->second->DisableHead(head);
    }
}
void GPU::bindVertexPuller(VertexPullerID vao)
{
    if (auto bufferit = InputAssemblers.find(vao); bufferit != InputAssemblers.end())
        activeIA = bufferit->second.get();
}
void GPU::bindVertexPuller(InputLayout* in)
{
    activeIA = in;
}
void GPU::unbindVertexPuller()
{
    activeIA = nullptr;
}
bool GPU::isVertexPuller(VertexPullerID vao)
{
    return InputAssemblers.find(vao) != InputAssemblers.end();
}

/// @}

ProgramID GPU::createProgram()
{
    auto ptr = std::make_unique<ShaderProgramm>();
    ProgramID id = reinterpret_cast<ProgramID>(ptr.get());
    ShaderStage.emplace(id, std::move(ptr));
    return id;
}
ProgramID GPU::createProgram(ShaderProgramm** _out, VertexShader vs, FragmentShader fs)
{
    auto ptr = std::make_unique<ShaderProgramm>(fs,vs);
    *_out = ptr.get();
    ProgramID id = reinterpret_cast<ProgramID>(ptr.get());
    ShaderStage.emplace(id, std::move(ptr));
    return id;
}
void GPU::deleteProgram(ProgramID prg)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
    {
        ShaderStage.erase(bufit);
    }
}
void GPU::attachShaders(ProgramID prg, VertexShader vs, FragmentShader fs)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
    {
        bufit->second->SetPixelShader(fs);
        bufit->second->SetVertexShader(vs);
    }
}
void GPU::setVS2FSType(ProgramID prg,uint32_t attrib,AttributeType type)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
    {
        bufit->second->SetVStoPSAttibutes(attrib, type);
    }
}
void GPU::useProgram(ProgramID prg)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
        activeSP = bufit->second.get();
}
void GPU::useProgram(ShaderProgramm* prg)
{
    activeSP = prg;
}
bool GPU::isProgram(ProgramID prg)
{
    return ShaderStage.find(prg) != ShaderStage.end();
}

void GPU::programUniform1f(ProgramID prg, uint32_t uniformId, float const&d)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
    {
        bufit->second->SetConstantBuffer(uniformId, d);
    }
}
void GPU::programUniform2f(ProgramID prg, uint32_t uniformId, glm::vec2 const&d)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
    {
        bufit->second->SetConstantBuffer(uniformId, d);
    }
}
void GPU::programUniform3f(ProgramID prg, uint32_t uniformId, glm::vec3 const&d)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
    {
        bufit->second->SetConstantBuffer(uniformId, d);
    }
}
void GPU::programUniform4f(ProgramID prg, uint32_t uniformId, glm::vec4 const&d)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
    {
        bufit->second->SetConstantBuffer(uniformId, d);
    }
}
void GPU::programUniformMatrix4f(ProgramID prg, uint32_t uniformId, glm::mat4 const&d)
{
    if (auto bufit = ShaderStage.find(prg); bufit != ShaderStage.end())
    {
        bufit->second->SetConstantBuffer(uniformId, d);
    }
}

/// @}

void GPU::createFramebuffer (uint32_t width_in, uint32_t height_in)
{
    width = width_in;
    height = height_in;
    size_t sz = (size_t)width * height;
    FrameBuffer.resize(sz);
    DepthStencil.resize(sz);
    const float HalfViewportWidth = width * 0.5f;
    const float HalfViewportHeight = height * 0.5f;
    Scale.v = DirectX::XMVectorSet(HalfViewportWidth, HalfViewportHeight, 1.0f, 0.0f);
    Offset.v = DirectX::XMVectorSet(HalfViewportWidth, HalfViewportHeight, 0.0f, 0.0f); //solve z!!!
    std::fill(DepthStencil.begin(), DepthStencil.end(), std::numeric_limits<float>::infinity());
}
void GPU::deleteFramebuffer()
{}
void GPU::resizeFramebuffer(uint32_t width_in, uint32_t height_in)
{
    createFramebuffer(width_in, height_in);
}
uint8_t* GPU::getFramebufferColor()
{
    return (uint8_t*)&FrameBuffer[0];
}
float* GPU::getFramebufferDepth()
{
    return &DepthStencil[0];
}
uint32_t GPU::getFramebufferWidth()
{
    return width;
}
uint32_t GPU::getFramebufferHeight()
{
    return height;
}

/// @}

void GPU::clear(float r,float g,float b,float a)
{
    namespace dxc = DirectX::PackedVector;
    dxc::XMCOLOR c; XMStoreColor2(&c, DirectX::XMVectorSet(r, g, b, a));
    std::fill(FrameBuffer.begin(), FrameBuffer.end(), c);
    std::fill(DepthStencil.begin(), DepthStencil.end(), std::numeric_limits<float>::infinity());
}


void GPU::drawTriangles(uint32_t nofVertices)
{
    if(activeIA->IsIndexed())
    {
        switch (activeIA->GetIndexType())
        {
        case IndexType::UINT8:
            drawIndexed<uint8_t>(nofVertices);
            return;
        case IndexType::UINT16:
            drawIndexed<uint16_t>(nofVertices);
            return;
        case IndexType::UINT32:
            drawIndexed<uint32_t>(nofVertices);
            return;
        default:
            return;
        }
    }

    auto& Assemler = *activeIA;
    auto& Program = *activeSP;
    std::vector<XMVSOut> VSOut{ size_t(nofVertices) };
    uint32_t SV_VertexID = 0;

    for(auto& Vertex : VSOut)
    {
        Program.InvokeVS(Vertex.gl, Assemler.MakeVertex(SV_VertexID));
        SV_VertexID++;
    }
    AssembleTriangles(VSOut);
}
void GPU::AssembleTriangles(std::vector<XMVSOut>& VSOut)
{
    for (size_t it = 0u; it < VSOut.size(); it += 3)
    {
        auto& v0 = VSOut[it];
        auto& v1 = VSOut[it + 1];
        auto& v2 = VSOut[it + 2];

        ClipCullTriangles(v0, v1, v2);
    }
}
void GPU::ClipCullTriangles(XMVSOut& v0, XMVSOut& v1, XMVSOut& v2)
{
    using namespace DirectX; 
    const auto Clip1V = [this](XMVSOut& v0, XMVSOut& v1, XMVSOut& v2)
    {
        auto vosize = activeSP->GetMonotonicSize();

        const float alphaA = (-v0.dx.SV_Position.f[3] - v0.dx.SV_Position.f[2]) / (v1.dx.SV_Position.f[3] - v0.dx.SV_Position.f[3] + v1.dx.SV_Position.f[2] - v0.dx.SV_Position.f[2]);
        const float alphaB = (-v0.dx.SV_Position.f[3] - v0.dx.SV_Position.f[2]) / (v2.dx.SV_Position.f[3] - v0.dx.SV_Position.f[3] + v2.dx.SV_Position.f[2] - v0.dx.SV_Position.f[2]);

        // interpolate to get v0a and v0b
        auto v0a = VSOutInterpolate(v0, v1, alphaA, vosize);
        auto v0b = VSOutInterpolate(v0, v2, alphaB, vosize);
        auto v0a2 = v0a;
        auto v2b = v2;
        // draw triangles
        PostProcessTriangle(v0a, v1, v2);
        PostProcessTriangle(v0b, v0a2, v2b);
    };
    const auto Clip2V = [this](XMVSOut& v0, XMVSOut& v1, XMVSOut& v2)
    {
        auto vosize = activeSP->GetMonotonicSize();
        // calculate alpha values for getting adjusted vertices
        const float alpha0 = (-v0.dx.SV_Position.f[3] - v0.dx.SV_Position.f[2]) / (v2.dx.SV_Position.f[3] - v0.dx.SV_Position.f[3] + v2.dx.SV_Position.f[2] - v0.dx.SV_Position.f[2]);
        const float alpha1 = (-v1.dx.SV_Position.f[3] - v1.dx.SV_Position.f[2]) / (v2.dx.SV_Position.f[3] - v1.dx.SV_Position.f[3] + v2.dx.SV_Position.f[2] - v1.dx.SV_Position.f[2]);
        // interpolate to get v0a and v0b

        v0 = VSOutInterpolate(v0, v2, alpha0, vosize);
        v1 = VSOutInterpolate(v1, v2, alpha1, vosize);
        // draw triangles
        PostProcessTriangle(v0, v1, v2);
    };

    VMFLOAT32A V0 = v0.dx.SV_Position;
    VMFLOAT32A V1 = v1.dx.SV_Position;
    VMFLOAT32A V2 = v2.dx.SV_Position;

    // Compare againgst W value
    XMVECTOR CT0 = XMVectorSplatW(V0);
    XMVECTOR CT1 = XMVectorSplatW(V1);
    XMVECTOR CT2 = XMVectorSplatW(V2);

    XMVECTOR R01 = XMVectorLess(CT0, V0);
    XMVECTOR R11 = XMVectorLess(CT1, V1);
    XMVECTOR R21 = XMVectorLess(CT2, V2);

    XMVECTOR RR1 = XMVectorAndInt(XMVectorAndInt(R01, R11), R21);
    if (_mm_movemask_ps(RR1) != 0) return;

    CT0 = XMVectorNegate(CT0);
    CT1 = XMVectorNegate(CT1);
    CT2 = XMVectorNegate(CT2);

    XMVECTOR R02 = XMVectorLess(V0, CT0);
    XMVECTOR R12 = XMVectorLess(V1, CT1);
    XMVECTOR R22 = XMVectorLess(V2, CT2);
    XMVECTOR RR2 = XMVectorAndInt(XMVectorAndInt(R02, R12), R22);
    if (_mm_movemask_ps(RR2) != 0) return;

    RR1 = XMVectorMergeZW(R02, R12);
    RR2 = _mm_shuffle_ps(RR1, R22, _MM_SHUFFLE(3, 2, 1, 0));

    int selector = _mm_movemask_ps(RR2) & 7;
    switch (selector)
    {
    case 0: PostProcessTriangle(v0, v1, v2); break;
    case 1: Clip1V(v0, v1, v2); break;
    case 2: Clip1V(v1, v2, v0); break;
    case 3: Clip2V(v0, v1, v2); break;
    case 4: Clip1V(v2, v0, v1); break;
    case 5: Clip2V(v2, v0, v1); break;
    case 6: Clip2V(v1, v2, v0); break;
    }
}
void GPU::PostProcessTriangle(XMVSOut& v0, XMVSOut& v1, XMVSOut& v2)
{
    using namespace DirectX;
    // homo -> NDC space transformation
    XMVECTOR wInv = XMVectorReciprocal(_mm_shuffle_ps(XMVectorMergeZW(v0.dx.SV_Position, v1.dx.SV_Position), v2.dx.SV_Position, _MM_SHUFFLE(3, 3, 3, 2)));
    XMMATRIX X;
    X.r[0] = XMVectorSplatX(wInv); // 1/w0
    X.r[1] = XMVectorSplatY(wInv); // 1/w1
    X.r[2] = XMVectorSplatZ(wInv); // 1/w2

    //Screen space transform, 1/w is stored in W
    v0.dx.SV_Position = XMVectorMultiplyAdd(XMVectorMultiply(v0.dx.SV_Position, X.r[0]), Scale, Offset); XMStoreFloat(&v0.dx.SV_Position.f[3], X.r[0]);
    v1.dx.SV_Position = XMVectorMultiplyAdd(XMVectorMultiply(v1.dx.SV_Position, X.r[1]), Scale, Offset); XMStoreFloat(&v1.dx.SV_Position.f[3], X.r[1]);
    v2.dx.SV_Position = XMVectorMultiplyAdd(XMVectorMultiply(v2.dx.SV_Position, X.r[2]), Scale, Offset); XMStoreFloat(&v2.dx.SV_Position.f[3], X.r[2]);
    
    //cull backfaces
    if (VMVector3Cross((v1.dx.SV_Position.v - v0.dx.SV_Position.v), (v2.dx.SV_Position.v - v0.dx.SV_Position.v)).f[2] < 0.0f)
        return;

    auto vosize = activeSP->GetMonotonicSize();
    for (unsigned i = 0; i < vosize; i++)
    {
        v0.dx.attributes[i].v = v0.dx.attributes[i].v * X.r[0];
        v1.dx.attributes[i].v = v1.dx.attributes[i].v * X.r[1];
        v2.dx.attributes[i].v = v2.dx.attributes[i].v * X.r[2];
    }
    DrawTriangle(v0, v1, v2);
}

void GPU::DrawTriangle(const XMVSOut& Vo0, const XMVSOut& Vo1, const XMVSOut& Vo2)
{
    // using pointers so we can swap (for sorting purposes)
    const auto* pv0 = &Vo0;
    const auto* pv1 = &Vo1;
    const auto* pv2 = &Vo2;

    // sorting vertices by y
    if (pv1->dx.SV_Position.f[1] < pv0->dx.SV_Position.f[1]) std::swap(pv0, pv1);
    if (pv2->dx.SV_Position.f[1] < pv1->dx.SV_Position.f[1]) std::swap(pv1, pv2);
    if (pv1->dx.SV_Position.f[1] < pv0->dx.SV_Position.f[1]) std::swap(pv0, pv1);

    if (pv0->dx.SV_Position.f[1] == pv1->dx.SV_Position.f[1]) // natural flat top
    {
        // sorting top vertices by x
        if (pv1->dx.SV_Position.f[0] < pv0->dx.SV_Position.f[0]) std::swap(pv0, pv1);

        DrawFlatTopTriangle(*pv0, *pv1, *pv2);
    }
    else if (pv1->dx.SV_Position.f[1] == pv2->dx.SV_Position.f[1]) // natural flat bottom
    {
        // sorting bottom vertices by x
        if (pv2->dx.SV_Position.f[0] < pv1->dx.SV_Position.f[0]) std::swap(pv1, pv2);

        DrawFlatBottomTriangle(*pv0, *pv1, *pv2);
    }
    else // general triangle
    {
        // find splitting vertex interpolant
        const float alphaSplit =
            (pv1->dx.SV_Position.f[1] - pv0->dx.SV_Position.f[1]) /
            (pv2->dx.SV_Position.f[1] - pv0->dx.SV_Position.f[1]);
        const auto vi = VSOutInterpolate(*pv0, *pv2, alphaSplit, activeSP->GetMonotonicSize());

        if (pv1->dx.SV_Position.f[0] < vi.dx.SV_Position.f[0]) // major right
        {
            DrawFlatBottomTriangle(*pv0, *pv1, vi);
            DrawFlatTopTriangle(*pv1, vi, *pv2);
        }
        else // major left
        {
            DrawFlatBottomTriangle(*pv0, vi, *pv1);
            DrawFlatTopTriangle(vi, *pv1, *pv2);
        }
    }
}
void GPU::DrawFlatTopTriangle(const XMVSOut& Vo0, const XMVSOut& Vo1, const XMVSOut& Vo2)
{
    using namespace DirectX;
    auto vosize = activeSP->GetMonotonicSize();
    XMVSOut InterpolantEdge = Vo1;
    // calulcate dVertex / dy
    // change in interpolant for every 1 change in y
    auto dv0 = XMVSOut::Subtract(Vo2, Vo0, vosize);
    auto dv1 = XMVSOut::Subtract(Vo2, Vo1, vosize);
    FXMVECTOR delta_Y = XMVectorSplatY(XMVectorReciprocal(dv0.dx.SV_Position));

    // delta over 0 and 1 resp
    dv0.Scale(delta_Y, vosize);
    dv1.Scale(delta_Y, vosize);

    // call the flat triangle render routine, right edge interpolant is it1
    DrawFlatTriangle(Vo0, Vo2, dv0, dv1, InterpolantEdge);
}
void GPU::DrawFlatBottomTriangle(const XMVSOut& Vo0, const XMVSOut& Vo1, const XMVSOut& Vo2)
{
    using namespace DirectX;
    auto vosize = activeSP->GetMonotonicSize();
    XMVSOut InterpolantEdge = Vo0;

    // calulcate dVertex / dy
    // change in interpolant for every 1 change in y
    auto dv0 = XMVSOut::Subtract(Vo1, Vo0, vosize);
    auto dv1 = XMVSOut::Subtract(Vo2, Vo0, vosize);
    FXMVECTOR delta_Y = XMVectorReciprocal(XMVectorSplatY(dv0.dx.SV_Position)); // minimize div (reciprocalEst is not good enough)

    // delta over 1 and 2 resp
    dv0.Scale(delta_Y, vosize);
    dv1.Scale(delta_Y, vosize);

    // call the flat triangle render routine, right edge interpolant is it0
    DrawFlatTriangle(Vo0, Vo2, dv0, dv1, InterpolantEdge);
}

void GPU::DrawFlatTriangle(const XMVSOut& it0,
    const XMVSOut& it2,
    const XMVSOut& dv0,
    const XMVSOut& dv1,
    XMVSOut& itEdge1)
{
    using namespace DirectX;
    const size_t Size = activeSP->GetMonotonicSize();

    // create edge interpolant for left edge (always v0)
    XMVSOut itEdge0 = it0;
    XMVSOut iLine;
    XMVSOut diLine;
    XMInPixel _P;

    // calculate start and end scanlines (AABB)
    const int yStart = std::max((int)ceil(it0.dx.SV_Position.f[1] - 0.5f), 0);
    const int yEnd = std::min((int)ceil(it2.dx.SV_Position.f[1] - 0.5f), (int)height); // the scanline AFTER the last line drawn

    // do interpolant prestep
    FXMVECTOR step = XMVectorReplicate(((float)yStart + 0.5f - it0.dx.SV_Position.f[1]));

    for (unsigned i = 0; i < Size; i++)
    {
        itEdge0.dx.attributes[i].v = XMVectorMultiplyAdd(dv0.dx.attributes[i], step, itEdge0.dx.attributes[i]);
        itEdge1.dx.attributes[i].v = XMVectorMultiplyAdd(dv1.dx.attributes[i], step, itEdge1.dx.attributes[i]);
    }
    itEdge0.dx.SV_Position.v = XMVectorMultiplyAdd(dv0.dx.SV_Position, step, itEdge0.dx.SV_Position);
    itEdge1.dx.SV_Position.v = XMVectorMultiplyAdd(dv1.dx.SV_Position, step, itEdge1.dx.SV_Position);

    for (int y = yStart; y < yEnd; y++, itEdge0.Increase(dv0, Size), itEdge1.Increase(dv1, Size))
    {
        // calculate start and end pixels
        const int xStart = std::max((int)ceil(itEdge0.dx.SV_Position.f[0] - 0.5f), 0);
        const int xEnd = std::min((int)ceil(itEdge1.dx.SV_Position.f[0] - 0.5f), (int)width); // the pixel AFTER the last pixel drawn

        // create scanline interpolant startpoint
        // (some waste for interpolating x,y,z, but makes life easier not having
        //  to split them off, and z will be needed in the future anyways...)

        iLine = itEdge0;
        FXMVECTOR step2 = XMVectorReplicate((float)xStart + 0.5f - itEdge0.dx.SV_Position.f[0]);
        FXMVECTOR Delta_X = XMVectorReciprocal(XMVectorSplatX(itEdge1.dx.SV_Position.v - itEdge0.dx.SV_Position.v));
        diLine = XMVSOut::Multiply(XMVSOut::Subtract(itEdge1, iLine, Size), Delta_X, Size);
        iLine += XMVSOut::Multiply(diLine, step2, Size);

        const size_t premulI = y * size_t(width);


        for (int x = xStart; x < xEnd; x++, iLine.Increase(diLine, Size))
        {
            if (auto [pass, zv] = DepthTest(x, premulI, iLine.dx.SV_Position.f[2]); pass)
            {
                // recover interpolated z from interpolated 1/z
                FXMVECTOR w = XMVectorReciprocalEst(XMVectorSplatW(iLine.dx.SV_Position));
                for (unsigned i = 0; i < Size; i++)
                    _P.dx.dx.attributes[i].v = iLine.dx.attributes[i].v * w;
                _P.dx.dx.SV_Position.v = iLine.dx.SV_Position;
                // invoke pixel shader with interpolated vertex attributes
                // and use result to set the pixel color on the screen
                DirectX::PackedVector::XMCOLOR col;
                XMStoreColor2(&col, activeSP->InvokePS(_P.gl).SV_Target);
                if (zv == DepthStencil[premulI + x]) //afxmt sanity check for data races
                {
                    FrameBuffer[premulI + x] = col;
                }

            }
        }
    }
}
std::pair<bool, float> GPU::DepthTest(uint32_t width_in, size_t PremulIndex, float z)
{
    //Warning! lockless programming!
    auto* reg = reinterpret_cast<std::atomic<float>*>(&DepthStencil[PremulIndex + width_in]);
    auto zv = std::atomic_load(reg);
    
    do
    {
        if (z >= zv) return { false, z };
    } while (!std::atomic_compare_exchange_weak(reg, &zv, z));

    return { true, z };
}

/// @}
