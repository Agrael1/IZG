#pragma once
#include <student/method.hpp>

class PhongMethod: public Method
{
public:
    PhongMethod();
    virtual ~PhongMethod();
public:
    virtual void onDraw(glm::mat4 const&proj,glm::mat4 const&view,glm::vec3 const&light,glm::vec3 const&camera) override;
private:
    ShaderProgramm* Programm;
    Buffer* VertexBuffer;
    Buffer* IndexBuffer;
    InputLayout* pInputLayout;
};
