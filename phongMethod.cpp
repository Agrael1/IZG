#include <student/phongMethod.hpp>
#include <student/bunny.hpp>


DirectX::XMVECTOR ProcTexWave(VMFLOAT32A wpos)
{
    wpos.f[0] += sinf((wpos.f[1] * 10)) / 10;
    wpos.f[0] *= 5;
    if (floorf(wpos.f[0] - floorf(wpos.f[0]) - 0.5f))
    {
        return DirectX::XMVectorSet(0.0f, 0.5f, 0.0f, 0.0f);
    }
    else
    {
        return  DirectX::XMVectorSet(1.0f, 1.0f, 0.0f, 0.0f);
    }
}
void phong_VS(OutVertex& outVertex, const InVertex& inVertex, const Uniforms& uniforms)
{
    using namespace DirectX;
    outVertex.gl_Position = outVertex.attributes[0].v4 = glm::vec4(inVertex.attributes[0].v3, 1.0f);
    outVertex.attributes[1] = inVertex.attributes[1];

    if (uniforms.uniform[4].v1 == std::numeric_limits<float>().infinity())
    {
        const auto mvpx = XMLoadFloat4x4((XMFLOAT4X4*)&uniforms.uniform[5]);
        XMStoreFloat4((XMFLOAT4*)&outVertex.gl_Position, XMVector4Transform(DirectX::XMLoadFloat4((XMFLOAT4*)&outVertex.gl_Position), mvpx));
        return;
    }

    const auto viewx = XMLoadFloat4x4((XMFLOAT4X4*)&uniforms.uniform[0].m4);
    const auto projx = XMLoadFloat4x4((XMFLOAT4X4*)&uniforms.uniform[1].m4);
    const auto mvpx = viewx*projx;
    XMStoreFloat4((XMFLOAT4*)&outVertex.gl_Position, XMVector4Transform(DirectX::XMLoadFloat4((XMFLOAT4*)&outVertex.gl_Position), mvpx));
}
void phong_FS(OutFragment& outFragment, const InFragment& inFragment, const Uniforms& uniforms)
{
    using namespace DirectX;
    const VMFLOAT32A wpos = XMLoadFloat4((XMFLOAT4*)&inFragment.attributes[0].v2.x);
    const VMFLOAT32A surf_norm = XMVector3Normalize(XMLoadFloat4((XMFLOAT4*)&inFragment.attributes[1].v4.x));

    const VMFLOAT32A lightDir = XMVector3Normalize(XMLoadFloat4((XMFLOAT4*)&uniforms.uniform[2].v4.x) - wpos.v);
    const VMFLOAT32A viewDir = XMVector3Normalize(XMLoadFloat4((XMFLOAT4*)&uniforms.uniform[3].v4.x) - wpos.v);

    auto objectColor = ProcTexWave(wpos);
    if(surf_norm.f[1] > 0.0f)
        objectColor = XMVectorLerp(objectColor, g_XMOne, surf_norm.f[1] * surf_norm.f[1]);

    auto diffuse = XMVectorMax(DirectX::XMVectorZero(), XMVector3Dot(surf_norm, lightDir));
    auto reflectDir = XMVector3Reflect(-lightDir.v, surf_norm);
    auto specular = XMVectorPow(XMVectorMax(DirectX::XMVectorZero(), XMVector3Dot(viewDir, reflectDir)), DirectX::XMVectorReplicate(40.0f));
    auto result = XMVectorSaturate(XMVectorMultiplyAdd(diffuse, objectColor, specular));

    XMStoreFloat4((XMFLOAT4*)&outFragment.gl_FragColor, result);
}


PhongMethod::PhongMethod()
{
    gpu.createProgram(&Programm, phong_VS, phong_FS);
    Programm->SetVStoPSAttibutes({ AttributeType::VEC3, AttributeType::VEC3 });

    constexpr auto const verticesSize = sizeof(BunnyVertex) * 1048;
    gpu.createBuffer(&VertexBuffer, verticesSize);
    VertexBuffer->InsertDataChunk(bunnyVertices, 0, verticesSize);

    gpu.createBuffer(&IndexBuffer, sizeof(bunnyIndices));
    IndexBuffer->InsertDataChunk(bunnyIndices, 0, sizeof(bunnyIndices));

    gpu.createVertexPuller(&pInputLayout);
    pInputLayout->SetDescriptor(0, AttributeType::VEC3, sizeof(BunnyVertex), 0, VertexBuffer);
    pInputLayout->SetDescriptor(1, AttributeType::VEC3, sizeof(BunnyVertex), sizeof(glm::vec3), VertexBuffer);
    pInputLayout->SetIndexing(IndexBuffer, IndexType::UINT32);
    pInputLayout->EnableHead(0);
    pInputLayout->EnableHead(1);
}
PhongMethod::~PhongMethod()
{
    gpu.deleteVertexPuller((VertexPullerID)pInputLayout);
    gpu.deleteBuffer(VertexBuffer->GetID());
    gpu.deleteBuffer(IndexBuffer->GetID());
    gpu.deleteProgram((ProgramID)Programm);
}

void PhongMethod::onDraw(const glm::mat4& proj, const glm::mat4& view, const glm::vec3& light, const glm::vec3& camera)
{
    gpu.clear(.5f, .5f, .5f, 1.f);

    gpu.bindVertexPuller(pInputLayout);
    gpu.useProgram(Programm);

    const glm::mat4& viewproj = proj * view;
    const float x(std::numeric_limits<float>().infinity());
    Programm->SetConstantBuffer(0, view, proj, light, camera, x, viewproj);

    gpu.drawIndexed(2092*3);
    gpu.unbindVertexPuller();
}

