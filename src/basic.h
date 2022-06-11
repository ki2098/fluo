#ifndef _BASIC_H
#define _BASIC_H

#define idx3d(i,j,k,size) ((i)*size[1]*size[2] + (j)*size[2]+(k))
#define idx4d(i,j,k,m,size) ((i)*size[1]*size[2]*size[3] + (j)*size[2]*size[3] + (k)*size[3] + (m))
#define idx5d(i,j,k,m,n,size) ((i)*size[1]*size[2]*size[3]*size[4] + (j)*size[2]*size[3]*size[4] + (k)*size[3]*size[4] + (m)*size[4] +(n))

static const int GUIDE = 2;

#endif