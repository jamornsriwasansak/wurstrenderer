#pragma once

#include <iostream>

#include <atomic>
#include <algorithm>
#include <cassert>
#include <exception>
#include <fstream>
#include <functional>
#include <future>
#include <iomanip>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#if __has_include(<filesystem>)
	#include <filesystem>
	namespace filesystem = std::filesystem;
#elif __has_include(<experimental/filesystem>)
	#include <experimental/filesystem>
	namespace filesystem = std::experimental::filesystem;
#else
	#include <ext/ghc/filesystem.hpp>
	namespace filesystem = ghc::filesystem;
#endif

#ifdef WURST_UNION_AS_STRUCT
#define SUnion struct
#else
#define SUnion union
#endif

#ifdef WURST_LIBTORCH
// libtorch stuffs
#pragma warning(push, 0) 
#include <torch/csrc/autograd/variable.h>
#include <torch/csrc/autograd/function.h>
#include <torch/types.h>
#include <torch/cuda.h>
#include <torch/data.h>
#include <torch/jit.h>
#include <torch/optim.h>
#include <torch/serialize.h>
#include <torch/types.h>
#include <torch/utils.h>
#pragma warning(pop)
#endif

// common stuffs
using std::make_pair;
using std::shared_ptr;
using std::unique_ptr;
using std::weak_ptr;
using std::const_pointer_cast;
using std::dynamic_pointer_cast;
//using std::reinterpret_pointer_cast;
using std::static_pointer_cast;
using std::make_shared;
using std::make_unique;

#include "common/math/spectrum.h"
using Spectrum = RgbSpectrum;

#include "common/math/math.h"

struct Bsdf;
struct Camera;
struct Camera2;
struct CdfTable;
struct Film;
struct Geometry;
struct Light;
struct Medium;
struct PhaseFunction;
struct ProgressReport;
struct Renderer;
struct Rng;
struct Sampler;
struct Scene;
struct Subpath;
struct SubpathSampler;
struct SubpathSampler2;
struct StopWatch;
struct Vertex;
struct Viewer;
struct Visualizer2;

template <typename> struct Fimage;
template <typename> struct SphRep;
template <typename> struct Texture;

struct LightImaginaryVertex;
struct LightVertex;
struct EnvmapVertex;
struct SurfaceVertex;
struct MediumVertex;
struct CameraVertex;
struct CameraImaginaryVertex;

#include "ext/json/json.hpp"
using nlohmann::json;
