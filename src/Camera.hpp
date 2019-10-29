/*
  Copyright 2015 - UVSQ
  Authors list: Nathalie Möller, Eric Petit

  This file is part of the FMM-viz.

  FMM-viz is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later version.

  FMM-viz is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License along with
  the FMM-viz. If not, see <http://www.gnu.org/licenses/>.
*/

#include <glm/glm.hpp>


class Camera
{
public:
    enum Direction { LEFT, RIGHT, UP, DOWN, FORWARD, BACK };

private:
  
  glm::vec3 position;
  glm::vec3 direction;
  glm::vec3 lookat;
  glm::vec3 up;

  float yaw_angle;
  float pitch_angle;
  
  float step;
  
  float fovy;
  
  glm::mat4 view;


public:
  Camera();

  void setPosition( glm::vec3 position );
  void setDirection( glm::vec3 position );

  void setStep( float step );
  
  void move( Direction direction );

  void yaw( float rad );

  void pitch( float rad );

  glm::vec3 getPosition();

  glm::vec3 getDirection();

  glm::vec3 getUp();

  void setMoveStep( float step );

private:
  void update();
  
};
