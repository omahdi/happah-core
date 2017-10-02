// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

namespace happah {

//DECLARATIONS

class RectangularCuboidLayout;

inline hpmat4x4 make_model_view_matrix(const RectangularCuboidLayout& layout, const hpmat4x4& viewMatrix, hpint i, hpint j, hpint k);

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

inline hpmat4x4 make_model_view_matrix(const RectangularCuboidLayout& layout, const hpmat4x4& viewMatrix, hpint i, hpint j, hpint k) {
     auto& lengths = layout.getLengths();
     auto& padding = layout.getPadding();
     auto x = i * (lengths.x + padding.x);
     auto y = j * (lengths.y + padding.y);
     auto z = k * (lengths.z + padding.z);
     return glm::translate(viewMatrix, Vector3D(x, y, z));
}

inline RectangularCuboidLayout make_rectangular_cuboid_layout(Point3D lengths, Point3D padding) { return { std::move(lengths), std::move(padding) }; }

}//namespace happah

