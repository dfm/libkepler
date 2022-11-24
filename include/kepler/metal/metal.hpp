#ifndef KEPLER_METAL_MATAL_HPP
#define KEPLER_METAL_MATAL_HPP

#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <string>

namespace kepler {
namespace metal {

// void solve();

enum SolverStatus {
  Uninitialized,
  OK,
  FailedToFindLibrary,
  FailedToFindFunction,
  FailedToCreatePipelineState,
  FailedToCreateCommandQueue,
  Success,
  FailedToCreateCommandBuffer,
  FailedToCreateCommandEncoder,
  UninitializedData,
};

template <typename T>
inline NS::String *function_name();

template <>
inline NS::String *function_name<float>() {
  return NS::String::string("solve_float", NS::ASCIIStringEncoding);
}

template <>
inline NS::String *function_name<double>() {
  return NS::String::string("solve_double", NS::ASCIIStringEncoding);
}

template <typename T>
struct solver {
  SolverStatus status = SolverStatus::Uninitialized;
  MTL::Device *device;
  MTL::ComputePipelineState *pipeline_state;
  MTL::CommandQueue *command_queue;

  bool initialized = false;
  std::size_t size = 0;
  MTL::Size grid_size, group_size;
  MTL::Buffer *eccentricity_buffer, *mean_anomaly_buffer, *eccentric_anomaly_buffer;

  void init(MTL::Device *device, MTL::Library *library) {
    this->device = device;

    if (!library) {
      status = SolverStatus::FailedToFindLibrary;
      return;
    }

    auto func_name = function_name<T>();
    auto function = library->newFunction(func_name);
    if (!function) {
      status = SolverStatus::FailedToFindFunction;
      return;
    }

    NS::Error *error;
    pipeline_state = this->device->newComputePipelineState(function, &error);
    if (!pipeline_state) {
      status = SolverStatus::FailedToCreatePipelineState;
      return;
    }

    command_queue = this->device->newCommandQueue();
    if (!command_queue) {
      status = SolverStatus::FailedToCreateCommandQueue;
      return;
    }

    status = SolverStatus::OK;
  }

  solver(MTL::Device *device, std::string library_path) {
    NS::Error *error;
    auto path = NS::String::string(library_path.c_str(), NS::StringEncoding::UTF8StringEncoding);
    auto library = device->newLibrary(path, &error);
    init(device, library);
  }

  solver(MTL::Device *device) {
    auto default_library = device->newDefaultLibrary();
    init(device, default_library);
  }

  void set_data(const T &eccentricity, std::size_t size, const T *mean_anomaly) {
    this->size = size;
    const auto buffer_size = size * sizeof(T);
    eccentricity_buffer = device->newBuffer(sizeof(T), MTL::ResourceStorageModeShared);
    mean_anomaly_buffer = device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);
    eccentric_anomaly_buffer = device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);
    ((T *)mean_anomaly_buffer->contents())[0] = eccentricity;

    auto target = (T *)mean_anomaly_buffer->contents();
    for (std::size_t i = 0; i < size; i++) {
      target[i] = mean_anomaly[i];
    }

    auto threads = pipeline_state->maxTotalThreadsPerThreadgroup();
    if (threads > size) {
      threads = size;
    }
    grid_size = MTL::Size(size, 1, 1);
    group_size = MTL::Size(threads, 1, 1);
    initialized = true;
  }

  SolverStatus execute() const {
    if (!initialized) {
      return SolverStatus::UninitializedData;
    }

    auto command_buffer = command_queue->commandBuffer();
    if (!command_buffer) {
      return SolverStatus::FailedToCreateCommandBuffer;
    }

    auto command_encoder = command_buffer->computeCommandEncoder();
    if (!command_encoder) {
      return SolverStatus::FailedToCreateCommandEncoder;
    }

    command_encoder->setComputePipelineState(pipeline_state);
    command_encoder->setBuffer(eccentricity_buffer, 0, 0);
    command_encoder->setBuffer(mean_anomaly_buffer, 0, 1);
    command_encoder->setBuffer(eccentric_anomaly_buffer, 0, 2);
    command_encoder->dispatchThreads(grid_size, group_size);
    command_encoder->endEncoding();

    command_buffer->commit();
    command_buffer->waitUntilCompleted();

    return SolverStatus::Success;
  }

  void get_results(T *eccentric_anomaly) const {
    auto source = (T *)eccentric_anomaly_buffer->contents();
    for (std::size_t i = 0; i < size; i++) {
      eccentric_anomaly[i] = source[i];
    }
  }

  SolverStatus solve(const T &eccentricity, std::size_t size, const T *mean_anomaly,
                     T *eccentric_anomaly) {
    set_data(eccentricity, size, mean_anomaly);
    auto status = execute();
    if (status == SolverStatus::Success) {
      get_results(eccentric_anomaly);
    }
    return status;
  }
};

}  // namespace metal
}  // namespace kepler

#endif