#include "Render.h"

#include "ospray/ospray_cpp.h"
#include "ospray/ospray_util.h"
#include "ospcommon/math/vec.h"
#include "Open3D/Open3D.h"
#include "Parameters.h"
#include <memory>
#include <vector>

using namespace ospcommon::math;

void Render::render(std::shared_ptr<open3d::geometry::TriangleMesh>& mesh, const Parameters& params, open3d::geometry::Image& imageOpen3D)
{
  vec2i imgSize;
  imgSize.x = imageOpen3D.width_; 
  imgSize.y = imageOpen3D.height_; 

  // Camera data
  Eigen::Vector3f epos(1.0f, -1.9f, 1.4f );
  vec3f cam_pos{ epos[0], epos[1], epos[2] };
  epos.normalize();
  Eigen::Vector3f dir = (-1.f) * epos;
  vec3f cam_view{ dir[0], dir[1], dir[2] };
  Eigen::Vector3f eright(dir[1], -dir[0], 0.f);
  Eigen::Vector3f eup = eright.cross(dir);
  vec3f cam_up{ eup[0], eup[1], eup[2] };

  // Triangle mesh data
  std::vector<vec3f> vertex;
  vertex.resize(mesh->vertices_.size());
  for (size_t i = 0; i < mesh->vertices_.size(); i++)
    vertex[i] = vec3f{ (float)mesh->vertices_[i][0], (float)mesh->vertices_[i][1], (float)mesh->vertices_[i][2] };
  std::vector<vec3f> normal;
  normal.resize(mesh->vertex_normals_.size());
  for (size_t i = 0; i < mesh->vertex_normals_.size(); i++)
    normal[i] = vec3f{ (float)mesh->vertex_normals_[i][0], (float)mesh->vertex_normals_[i][1], (float)mesh->vertex_normals_[i][2] };
  std::vector<vec3ui> index;
  index.resize(mesh->triangles_.size());
  for (size_t i = 0; i < mesh->triangles_.size(); i++)
    index[i] = vec3ui{ (unsigned int)mesh->triangles_[i][0], (unsigned int)mesh->triangles_[i][1], (unsigned int)mesh->triangles_[i][2] };

  // Colors data
  Eigen::Vector3f diffuseColor{ 0.95f, 0.9f, 0.9f };
  Eigen::Vector3f lightColor{ 1.0f, 1.0f, 1.0f };
  Eigen::Vector3f backgroundColor{ 0.85f, 0.9f, 0.95f };

  // Render with OSPRAY
  OSPError init_error = ospInit();
  if (init_error != OSP_NO_ERROR)
  {
    std::cout << "Ospray init error" << init_error << std::endl;
    return;
  }

  {
    // create and setup camera
    ospray::cpp::Camera camera("perspective");
    camera.setParam("aspect", imgSize.x / (float)imgSize.y);
    camera.setParam("position", cam_pos);
    camera.setParam("direction", cam_view);
    camera.setParam("up", cam_up);
    camera.commit(); // commit each object to indicate modifications are done

    ospray::cpp::World world;
    
    if (params.renderingQuality == RenderingQuality::Low)
    {
      // SCIVIS
      if (!vertex.empty())
      {
        // create and setup model and mesh
        ospray::cpp::Geometry mesh("mesh");
        mesh.setParam("vertex.position", ospray::cpp::Data(vertex));
        mesh.setParam("vertex.normal", ospray::cpp::Data(normal));
        mesh.setParam("index", ospray::cpp::Data(index));
        mesh.commit();

        ospray::cpp::Material material = ospNewMaterial("scivis", "obj");
        material.setParam("kd", vec3f(diffuseColor[0], diffuseColor[1], diffuseColor[2]));
        //material.setParam("ks", vec3f(0.7f, 0.7f, 0.6f));
        //material.setParam("ns", 25.f);
        material.setParam("d", 1.f);
        material.setParam("tf", 0.f);
        material.commit();

        // put the mesh into a model
        ospray::cpp::GeometricModel model(mesh);
        ospSetObject(model.handle(), "material", material.handle());
        model.commit();

        // put the model into a group (collection of models)
        ospray::cpp::Group group;
        group.setParam("geometry", ospray::cpp::Data(model));
        group.commit();

        // put the group into an instance (give the group a world transform)
        ospray::cpp::Instance instance(group);
        instance.commit();

        // put the instance in the world
        world.setParam("instance", ospray::cpp::Data(instance));
      }

      // create and setup light for Ambient Occlusion
      ospray::cpp::Light ambient_light("ambient");
      ambient_light.setParam("color", vec3f(lightColor[0], lightColor[1], lightColor[2]));
      ambient_light.setParam("intensity", 1.0f);
      ambient_light.commit();
      world.setParam("light", ospray::cpp::Data(ambient_light));
      world.commit();

      // create renderer, choose Scientific Visualization renderer
      ospray::cpp::Renderer renderer("scivis");
      renderer.setParam("aoSamples", 32);
      renderer.setParam("pixelSamples", 8);
      renderer.setParam("backgroundColor", vec3f(backgroundColor[0], backgroundColor[1], backgroundColor[2]));
      renderer.commit();

      // create and setup framebuffer
      ospray::cpp::FrameBuffer framebuffer(
        imgSize, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
      framebuffer.clear();
      framebuffer.renderFrame(renderer, camera, world);

      uint32_t* fb = (uint32_t*)framebuffer.map(OSP_FB_COLOR);

      // convert resulting frame buffer in an Open3D image
      int k = 0;
      for (int i = 0; i < imageOpen3D.height_; i++)
        for (int j = 0; j < imageOpen3D.width_; j++)
        {
          unsigned char* tab = (unsigned char*)&fb[(imageOpen3D.height_ - 1 - i) * imageOpen3D.width_ + j];
          imageOpen3D.data_[k++] = tab[0];
          imageOpen3D.data_[k++] = tab[1];
          imageOpen3D.data_[k++] = tab[2];
        }
      framebuffer.unmap(fb);
    }
    else
    {
      // PATHTRACER
      if (!vertex.empty())
      {
        // create and setup model and mesh
        ospray::cpp::Geometry mesh("mesh");
        mesh.setParam("vertex.position", ospray::cpp::Data(vertex));
        mesh.setParam("vertex.normal", ospray::cpp::Data(normal));
        mesh.setParam("index", ospray::cpp::Data(index));
        mesh.commit();

        ospray::cpp::Material material = ospNewMaterial("pathtracer", "obj");
        material.setParam("kd", vec3f(diffuseColor[0], diffuseColor[1], diffuseColor[2]));
        //material.setParam("ks", vec3f(0.7f, 0.7f, 0.6f));
        //material.setParam("ns", 25.f);
        material.setParam("d", 1.f);
        material.setParam("tf", 0.f);
        material.commit();

        // put the mesh into a model
        ospray::cpp::GeometricModel model(mesh);
        ospSetObject(model.handle(), "material", material.handle());
        model.commit();

        // put the model into a group (collection of models)
        ospray::cpp::Group group;
        group.setParam("geometry", ospray::cpp::Data(model));
        group.commit();

        // put the group into an instance (give the group a world transform)
        ospray::cpp::Instance instance(group);
        instance.commit();

        // put the instance in the world
        world.setParam("instance", ospray::cpp::Data(instance));
      }

      ospray::cpp::Light sun_light("sunSky");
      sun_light.setParam("up", vec3f(0.f, 0.f, 1.f));
      sun_light.setParam("direction", vec3f(-0.2f, 0.2f, -0.8f));
      //sun_light.setParam("color", vec3f(lightColor[0], lightColor[1], lightColor[2]));
      //sun_light.setParam("intensity", 0.9f);
      sun_light.commit();
      world.setParam("light", ospray::cpp::Data(sun_light));

      world.commit();

      // create renderer
      ospray::cpp::Renderer renderer("pathtracer");
      renderer.setParam("roulettePathLength", 3);
      renderer.setParam("pixelSamples", 64);
      //renderer.setParam("backgroundColor", vec3f(backgroundColor[0], backgroundColor[1], backgroundColor[2]));
      renderer.commit();

      // create and setup framebuffer
      ospray::cpp::FrameBuffer framebuffer(
        imgSize, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
      framebuffer.clear();
      framebuffer.renderFrame(renderer, camera, world);

      uint32_t* fb = (uint32_t*)framebuffer.map(OSP_FB_COLOR);

      // convert resulting frame buffer in an Open3D image
      int k = 0;
      for (int i = 0; i < imageOpen3D.height_; i++)
        for (int j = 0; j < imageOpen3D.width_; j++)
        {
          unsigned char* tab = (unsigned char*)&fb[(imageOpen3D.height_ - 1 - i) * imageOpen3D.width_ + j];
          imageOpen3D.data_[k++] = tab[0];
          imageOpen3D.data_[k++] = tab[1];
          imageOpen3D.data_[k++] = tab[2];
        }
      framebuffer.unmap(fb);
    }
  }
  ospShutdown();
  
}

