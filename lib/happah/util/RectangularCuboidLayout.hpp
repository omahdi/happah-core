// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

namespace happah {

//DECLARATIONS

class RectangularCuboidLayout;

inline std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const RectangularCuboidLayout& layout, glm::ivec3 min, glm::ivec3 max);

inline hpmat4x4 make_model_matrix(const RectangularCuboidLayout& layout, const hpmat4x4& matrix, hpint i, hpint j, hpint k);

inline hpmat4x4 make_model_matrix(const RectangularCuboidLayout& layout, hpint i, hpint j, hpint k);

inline RectangularCuboidLayout make_rectangular_cuboid_layout(Point3D lengths, Point3D padding = Point3D(0));

//DEFINITIONS

class RectangularCuboidLayout {
public:
     RectangularCuboidLayout(Point3D lengths, Point3D padding)
          : m_lengths(std::move(lengths)), m_padding(std::move(padding)) {}

     auto& getLengths() const { return m_lengths; }

     auto& getPadding() const { return m_padding; }

private:
     Point3D m_lengths;
     Point3D m_padding;

};//RectangularCuboidLayout

inline std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const RectangularCuboidLayout& layout, glm::ivec3 min, glm::ivec3 max) {
     auto& lengths = layout.getLengths();
     auto& padding = layout.getPadding();
     auto x0 = min.x * (lengths.x + padding.x);
     auto y0 = min.y * (lengths.y + padding.y);
     auto z0 = min.z * (lengths.z + padding.z);
     auto x1 = max.x * (lengths.x + padding.x) + lengths.x;
     auto y1 = max.y * (lengths.y + padding.y) + lengths.y;
     auto z1 = max.z * (lengths.z + padding.z) + lengths.z;

     return std::make_tuple(Point3D(x0, y0, z0), Point3D(x1, y1, z1));
}

inline hpmat4x4 make_model_matrix(const RectangularCuboidLayout& layout, const hpmat4x4& matrix, hpint i, hpint j, hpint k) {
     auto& lengths = layout.getLengths();
     auto& padding = layout.getPadding();
     auto x = i * (lengths.x + padding.x);
     auto y = j * (lengths.y + padding.y);
     auto z = k * (lengths.z + padding.z);

     return glm::translate(matrix, Vector3D(x, y, z));
}

inline hpmat4x4 make_model_matrix(const RectangularCuboidLayout& layout, hpint i, hpint j, hpint k) { return make_model_matrix(layout, hpmat4x4(1), i, j, k); }

inline RectangularCuboidLayout make_rectangular_cuboid_layout(Point3D lengths, Point3D padding) { return { std::move(lengths), std::move(padding) }; }

}//namespace happah

