/*! @file GraphicsDebugData.h
 *  @brief Data for sim to graph
 */

#ifndef VISUALIZATION_DATA_H
#define VISUALIZATION_DATA_H

#define VISUALIZATION_MAX_PATH_POINTS 2000
#define VISUALIZATION_MAX_PATHS 10
#define VISUALIZATION_MAX_ITEMS 100
#define VISUALIZATION_MAX_CHEETAHS 0

#include "cppTypes.h"

struct SphereVisualization
{
	Vec3<float> position;
	Vec4<float> color;
	double radius;
};

struct BlockVisualization
{
	Vec3<float> dimension;
	Vec3<float> corner_position;
	Vec3<float> rpy;
	Vec4<float> color;
};

struct ArrowVisualization
{
	Vec3<float> base_position;
	Vec3<float> direction;
	Vec4<float> color;
	float head_width;
	float head_length;
	float shaft_width;
};

struct CheetahVisualization
{
	Vec12<float> q;
	Quat<float> quat;
	Vec3<float> p;
	Vec4<float> color;
};

struct PathVisualization
{
	size_t num_points = 0;
	Vec4<float> color;
	Vec3<float> position[VISUALIZATION_MAX_PATH_POINTS];
};

struct ConeVisualization
{
	Vec3<float> point_position;
	Vec3<float> direction;
	Vec4<float> color;
	double radius;
};


struct VisualizationData
{
	size_t num_paths = 0, num_arrows = 0, num_cones = 0, num_spheres = 0, num_blocks = 0;
	SphereVisualization spheres[VISUALIZATION_MAX_ITEMS];
	BlockVisualization blocks[VISUALIZATION_MAX_ITEMS];
	ArrowVisualization arrows[VISUALIZATION_MAX_ITEMS];
	ConeVisualization cones[VISUALIZATION_MAX_ITEMS];
	PathVisualization paths[VISUALIZATION_MAX_PATHS];
};

#endif