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

#include "Camera.hpp"
#define GLM_FORCE_RADIANS
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <cmath>
using namespace std;


Camera::Camera()
{
  position  = glm::vec3( 0.0f, 0.0f, 0.0f );
  direction = glm::vec3( 0.0f, 0.0f, 1.0f );
  up        = glm::vec3( 0.0f, 1.0f, 0.0f );
  fovy = 45.0f;
  step = 0.1f;
}


void Camera::setPosition( glm::vec3 position )
{
  this->position = position;
}

void Camera::setDirection( glm::vec3 direction )
{
  this->direction = direction;
}

void Camera::update()
{
  //direction = glm::normalize( lookat - position );
  //view = glm::lookAt( position, lookat, up );
}


void Camera::move( Direction dir )
{
  switch( dir ) {
  case UP:
    position += up * step;
    break;
  case DOWN:
    position -= up * step;
    break;
  case LEFT:
    position += glm::cross( direction, up ) * step;
    break;
  case RIGHT:
    position -= glm::cross( direction, up ) * step;
    break;
  case FORWARD:
    position -= direction * step;
    break;
  case BACK:
    position += direction * step;
    break;
  }
  update();
}


void Camera::yaw( float angle )
{
  this->yaw_angle = angle;
  glm::vec3 axis = glm::cross( direction, up );
  glm::quat pitch_quat = glm::angleAxis( pitch_angle, axis );
  glm::quat heading_quat = glm::angleAxis(yaw_angle, up);
  glm::quat temp = glm::cross(pitch_quat, heading_quat);
  temp = glm::normalize(temp);
  direction = glm::rotate(temp, direction);
}


void Camera::pitch( float angle )
{
  this->pitch_angle = angle;
  glm::vec3 axis = glm::cross( direction, up );
  glm::quat pitch_quat = glm::angleAxis( pitch_angle, axis );
  glm::quat heading_quat = glm::angleAxis(yaw_angle, up);
  glm::quat temp = glm::cross(pitch_quat, heading_quat);
  temp = glm::normalize(temp);
  direction = glm::rotate(temp, direction);
}


glm::vec3 Camera::getPosition()
{
  return position;
}

glm::vec3 Camera::getDirection()
{
  return direction;
}

glm::vec3 Camera::getUp()
{
  return up;
}

void Camera::setStep( float step )
{
  this->step = step;
}
