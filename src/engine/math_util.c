#include <ultra64.h>

#include "sm64.h"
#include "engine/graph_node.h"
#include "math_util.h"
#include "surface_collision.h"

#include "trig_tables.inc.c"
#include "trig_tables_fixed.inc.c"

// Variables for a spline curve animation (used for the flight path in the grand star cutscene)
Vec4s *gSplineKeyframe;
float gSplineKeyframeFraction;
int gSplineState;

float rsqrtf(float number);

// These functions have bogus return values.
// Disable the compiler warning.
#pragma GCC diagnostic push

#ifdef __GNUC__
#if defined(__clang__)
  #pragma GCC diagnostic ignored "-Wreturn-stack-address"
#else
  #pragma GCC diagnostic ignored "-Wreturn-local-addr"
#endif
#endif

/// Copy vector 'src' to 'dest'
void vec3f_copy(Vec3f dest, Vec3f src) {
	for(int i = 0; i < 3; i++) {
    	dest[i] = src[i];
	}
}
void vec3q_copy(Vec3q destq, Vec3q srcq) {
	for(int i = 0; i < 3; i++) {
    	destq[i] = srcq[i];
	}
}

/// Set vector 'dest' to (x, y, z)
ALWAYS_INLINE void vec3f_set(Vec3f dest, f32 x, f32 y, f32 z) {
    dest[0] = x;
    dest[1] = y;
    dest[2] = z;
}
ALWAYS_INLINE void vec3q_set(Vec3q destq, q32 xq, q32 yq, q32 zq) {
    destq[0] = xq;
    destq[1] = yq;
    destq[2] = zq;
}

/// Add vector 'a' to 'dest'
void vec3f_add(Vec3f dest, Vec3f a) {
	for(int i = 0; i < 3; i++) {
    	dest[i] += a[i];
	}
}
void vec3q_add(Vec3q destq, Vec3q aq) {
	for(int i = 0; i < 3; i++) {
    	destq[i] += aq[i];
	}
}

/// Make 'dest' the sum of vectors a and b.
void vec3f_sum(Vec3f dest, Vec3f a, Vec3f b) {
	for(int i = 0; i < 3; i++) {
    	dest[i] = a[i] + b[i];
	}
}

/// Copy vector src to dest
void vec3s_copy(Vec3s dest, Vec3s src) {
	for(int i = 0; i < 3; i++) {
    	dest[i] = src[i];
	}
}

/// Set vector 'dest' to (x, y, z)
ALWAYS_INLINE void vec3s_set(Vec3s dest, s16 x, s16 y, s16 z) {
	dest[0] = x;
	dest[1] = y;
	dest[2] = z;
}

/// Add vector a to 'dest'
void vec3s_add(Vec3s dest, Vec3s a) {
	for(int i = 0; i < 3; i++) {
    	dest[i] += a[i];
	}
}

/// Make 'dest' the sum of vectors a and b.
void vec3s_sum(Vec3s dest, Vec3s a, Vec3s b) {
	for(int i = 0; i < 3; i++) {
    	dest[i] = a[i] + b[i];
	}
}

/// Subtract vector a from 'dest'
void vec3s_sub(Vec3s dest, Vec3s a) {
	for(int i = 0; i < 3; i++) {
    	dest[i] -= a[i];
	}
}
void vec3q_sub(Vec3q destq, Vec3q aq) {
	for(int i = 0; i < 3; i++) {
    	destq[i] -= aq[i];
	}
}

/// Convert short vector a to float vector 'dest'
void vec3s_to_vec3f(Vec3f dest, Vec3s a) {
	for(int i = 0; i < 3; i++) {
    	dest[i] = a[i];
	}
}

/**
 * Convert float vector a to a short vector 'dest' by rounding the components
 * to the nearest integer.
 */
void vec3q_to_vec3s(Vec3s dest, Vec3q aq) {
	for(int i = 0; i < 3; i++) {
    	// add/subtract 0.5 in order to round to the nearest s32 instead of truncating
		q32 vq = aq[i];
		//if(vq < 0) {
		//	vq -= q(0.5);
		//} else {
		//	vq += q(0.5);
		//}
    	dest[i] = qtrunc(vq);
	}
}
void vec3q_to_vec3f(Vec3f dest, Vec3q aq) {
	for(int i = 0; i < 3; i++) {
    	dest[i] = qtof(aq[i]);
	}
}

void vec3f_to_vec3s(Vec3s dest, Vec3f a) {
	for(int i = 0; i < 3; i++) {
		// add/subtract 0.5 in order to round to the nearest s32 instead of truncating
		q32 vq = q(a[i]);
		if(vq < 0) {
			vq -= q(0.5);
		} else {
			vq += q(0.5);
		}
    	dest[i] = qtrunc(vq);
    }
}

void vec3f_to_vec3q(Vec3q destq, Vec3f a) {
	for(int i = 0; i < 3; i++) {
    	destq[i] = q(a[i]);
    }
}

/**
 * Set 'dest' the normal vector of a triangle with vertices a, b and c.
 * It is similar to vec3f_cross, but it calculates the vectors (c-b) and (b-a)
 * at the same time.
 */
void find_vector_perpendicular_to_planeq(Vec3q dest, Vec3q a, Vec3q b, Vec3q c) {
    dest[0] = qmul(b[1] - a[1], c[2] - b[2]) - qmul(c[1] - b[1], b[2] - a[2]);
    dest[1] = qmul(b[2] - a[2], c[0] - b[0]) - qmul(c[2] - b[2], b[0] - a[0]);
    dest[2] = qmul(b[0] - a[0], c[1] - b[1]) - qmul(c[0] - b[0], b[1] - a[1]);
}

/// Make vector 'dest' the cross product of vectors a and b.
void vec3f_cross(Vec3f dest, Vec3f a, Vec3f b) {
    dest[0] = a[1] * b[2] - b[1] * a[2];
    dest[1] = a[2] * b[0] - b[2] * a[0];
    dest[2] = a[0] * b[1] - b[0] * a[1];
}

/// Make vector 'dest' the cross product of vectors a and b.
void vec3q_cross(Vec3q destq, Vec3q aq, Vec3q bq) {
    destq[0] = qmul(aq[1], bq[2]) - qmul(bq[1], aq[2]);
    destq[1] = qmul(aq[2], bq[0]) - qmul(bq[2], aq[0]);
    destq[2] = qmul(aq[0], bq[1]) - qmul(bq[0], aq[1]);
}

/// Scale vector 'dest' so it has length 1
void vec3f_normalize(Vec3f dest) {
    //! Possible division by zero
    f32 invsqrt = rsqrtf(dest[0] * dest[0] + dest[1] * dest[1] + dest[2] * dest[2]);

    dest[0] *= invsqrt;
    dest[1] *= invsqrt;
    dest[2] *= invsqrt;
}
/// Scale vector 'dest' so it has length 1
void vec3q_normalize(Vec3q destq) {
    q32 invsqrtq = rsqrtq(qmul(destq[0], destq[0]) + qmul(destq[1], destq[1]) + qmul(destq[2], destq[2]));

    destq[0] = qmul(destq[0], invsqrtq);
    destq[1] = qmul(destq[1], invsqrtq);
    destq[2] = qmul(destq[2], invsqrtq);
}

#pragma GCC diagnostic pop

/**
 * Set mtx to a look-at matrix for the camera. The resulting transformation
 * transforms the world as if there exists a camera at position 'from' pointed
 * at the position 'to'. The up-vector is assumed to be (0, 1, 0), but the 'roll'
 * angle allows a bank rotation of the camera.
 */
void mtx_lookat(ShortMatrix* mtx, Vec3q fromq, Vec3q toq, s16 roll) {
    s32 dxi = qtrunc(toq[0] - fromq[0]);
    s32 dyi = qtrunc(toq[1] - fromq[1]);
    s32 dzi = qtrunc(toq[2] - fromq[2]);

    register s32 lengthi = -sqrtu(dxi * dxi + dzi * dzi);
    q32 dxq = lengthi? dxi * ONE / lengthi: 0;
    q32 dzq = lengthi? dzi * ONE / lengthi: 0;

    q32 yColYq = cosqs(roll);
    q32 xColYq = qmul(sinqs(roll), dzq);
    q32 zColYq = qmul(-sinqs(roll), dxq);

    lengthi = -sqrtu(dxi * dxi + dyi * dyi + dzi * dzi);
    q32 xColZq = lengthi? dxi * ONE / lengthi: 0;
    q32 yColZq = lengthi? dyi * ONE / lengthi: 0;
    q32 zColZq = lengthi? dzi * ONE / lengthi: 0;

    q32 xColXq = qmul(yColYq, zColZq) - qmul(zColYq, yColZq);
    q32 yColXq = qmul(zColYq, xColZq) - qmul(xColYq, zColZq);
    q32 zColXq = qmul(xColYq, yColZq) - qmul(yColYq, xColZq);

    q32 lengthq = sqrtq(qmul(xColXq, xColXq) + qmul(yColXq, yColXq) + qmul(zColXq, zColXq));
    if(lengthq) {
        xColXq = qdiv(xColXq, lengthq);
        yColXq = qdiv(yColXq, lengthq);
        zColXq = qdiv(zColXq, lengthq);
    }

    xColYq = qmul(yColZq, zColXq) - qmul(zColZq, yColXq);
    yColYq = qmul(zColZq, xColXq) - qmul(xColZq, zColXq);
    zColYq = qmul(xColZq, yColXq) - qmul(yColZq, xColXq);

    lengthq = sqrtq(qmul(xColYq, xColYq) + qmul(yColYq, yColYq) + qmul(zColYq, zColYq));
    if(lengthq) {
        xColYq = qdiv(xColYq, lengthq);
        yColYq = qdiv(yColYq, lengthq);
        zColYq = qdiv(zColYq, lengthq);
    }

    mtx->m[0][0] = xColXq;
    mtx->m[1][0] = yColXq;
    mtx->m[2][0] = zColXq;
    mtx->t[0] = qtrunc(-(qmul(fromq[0], xColXq) + qmul(fromq[1], yColXq) + qmul(fromq[2], zColXq)));
    mtx->m[0][1] = -xColYq;
    mtx->m[1][1] = -yColYq;
    mtx->m[2][1] = -zColYq;
    mtx->t[1] = -qtrunc(-(qmul(fromq[0], xColYq) + qmul(fromq[1], yColYq) + qmul(fromq[2], zColYq)));
    mtx->m[0][2] = -xColZq;
    mtx->m[1][2] = -yColZq;
    mtx->m[2][2] = -zColZq;
    mtx->t[2] = -qtrunc(-(qmul(fromq[0], xColZq) + qmul(fromq[1], yColZq) + qmul(fromq[2], zColZq)));
}

#include <port/gfx/gfx.h>

/**
 * Set 'dest' to a transformation matrix that turns an object to face the camera.
 * 'mtx' is the look-at matrix from the camera
 * 'position' is the position of the object in the world
 * 'angle' rotates the object while still facing the camera.
 */
void mtx_billboard(ShortMatrix* dest, ShortMatrix* mtx, Vec3s position, s16 angle) {
    dest->m[0][0] = cosqs(angle);
    dest->m[1][0] = sinqs(angle);
    dest->m[2][0] = 0;

    dest->m[0][1] = dest->m[1][0];
    dest->m[1][1] = -dest->m[0][0];
    dest->m[2][1] = 0;

    dest->m[0][2] = 0;
    dest->m[1][2] = 0;
    dest->m[2][2] = 1;

    dest->t[0] = qtrunc(mtx->m[0][0] * position[0] + mtx->m[0][1] * position[1] + mtx->m[0][2] * position[2]) + mtx->t[0];
    dest->t[1] = qtrunc(mtx->m[1][0] * position[0] + mtx->m[1][1] * position[1] + mtx->m[1][2] * position[2]) + mtx->t[1];
    dest->t[2] = qtrunc(mtx->m[2][0] * position[0] + mtx->m[2][1] * position[1] + mtx->m[2][2] * position[2]) + mtx->t[2];
}

/**
 * Set 'dest' to a transformation matrix that aligns an object with the terrain
 * based on the normal. Used for enemies.
 * 'upDir' is the terrain normal
 * 'yaw' is the angle which it should face
 * 'pos' is the object's position in the world
 */
void mtx_align_terrain_normal(ShortMatrix* destq, Vec3q upDirq, Vec3q posq, s16 yaw) {
    Vec3q lateralDirq;
    Vec3q leftDirq;
    Vec3q forwardDirq;

    vec3q_set(lateralDirq, sinqs(yaw), 0, cosqs(yaw));
    vec3q_normalize(upDirq);

    vec3q_cross(leftDirq, upDirq, lateralDirq);
    vec3q_normalize(leftDirq);

    vec3q_cross(forwardDirq, leftDirq, upDirq);
    vec3q_normalize(forwardDirq);

    destq->m[0][0] = leftDirq[0];
    destq->m[0][1] = leftDirq[1];
    destq->m[0][2] = leftDirq[2];

    destq->m[1][0] = upDirq[0];
    destq->m[1][1] = upDirq[1];
    destq->m[1][2] = upDirq[2];

    destq->m[2][0] = forwardDirq[0];
    destq->m[2][1] = forwardDirq[1];
    destq->m[2][2] = forwardDirq[2];

    destq->t[0] = qtrunc(posq[0]);
    destq->t[1] = qtrunc(posq[1]);
    destq->t[2] = qtrunc(posq[2]);
}

/**
 * Set 'mtx' to a transformation matrix that aligns an object with the terrain
 * based on 3 height samples in an equilateral triangle around the object.
 * Used for Mario when crawling or sliding.
 * 'yaw' is the angle which it should face
 * 'pos' is the object's position in the world
 * 'radius' is the distance from each triangle vertex to the center
 */
void mtx_align_terrain_triangle(ShortMatrix* mtxq, Vec3f posf, s16 yaw, s32 radius) {
	Vec3q posq;
	vec3f_to_vec3q(posq, posf);
    struct Surface *sp74;
    Vec3q point0q;
    Vec3q point1q;
    Vec3q point2q;
    Vec3q forwardq;
    Vec3q xColumnq;
    Vec3q yColumnq;
    Vec3q zColumnq;
    q32 avgYq;
    q32 minYq = -q(radius) * 3;

    point0q[0] = posq[0] + radius * sinqs(yaw + 0x2AAA) / ONE;
    point0q[2] = posq[2] + radius * cosqs(yaw + 0x2AAA) / ONE;
    point1q[0] = posq[0] + radius * sinqs(yaw + 0x8000) / ONE;
    point1q[2] = posq[2] + radius * cosqs(yaw + 0x8000) / ONE;
    point2q[0] = posq[0] + radius * sinqs(yaw + 0xD555) / ONE;
    point2q[2] = posq[2] + radius * cosqs(yaw + 0xD555) / ONE;

    point0q[1] = find_floorq(point0q[0], posq[1] + q(150), point0q[2], &sp74);
    point1q[1] = find_floorq(point1q[0], posq[1] + q(150), point1q[2], &sp74);
    point2q[1] = find_floorq(point2q[0], posq[1] + q(150), point2q[2], &sp74);

    if (point0q[1] - posq[1] < minYq) {
        point0q[1] = posq[1];
    }

    if (point1q[1] - posq[1] < minYq) {
        point1q[1] = posq[1];
    }

    if (point2q[1] - posq[1] < minYq) {
        point2q[1] = posq[1];
    }

    avgYq = (point0q[1] + point1q[1] + point2q[1]) / 3;

    vec3q_set(forwardq, sinqs(yaw), 0, cosqs(yaw));
    find_vector_perpendicular_to_planeq(yColumnq, point0q, point1q, point2q);
    vec3q_normalize(yColumnq);
    vec3q_cross(xColumnq, yColumnq, forwardq);
    vec3q_normalize(xColumnq);
    vec3q_cross(zColumnq, xColumnq, yColumnq);
    vec3q_normalize(zColumnq);

    mtxq->m[0][0] = xColumnq[0];
    mtxq->m[0][1] = xColumnq[1];
    mtxq->m[0][2] = xColumnq[2];
    mtxq->t[0] = qtrunc(posq[0]);

    mtxq->m[1][0] = yColumnq[0];
    mtxq->m[1][1] = yColumnq[1];
    mtxq->m[1][2] = yColumnq[2];
    mtxq->t[1] = qtrunc((avgYq < posq[1]) ? posq[1] : avgYq);

    mtxq->m[2][0] = zColumnq[0];
    mtxq->m[2][1] = zColumnq[1];
    mtxq->m[2][2] = zColumnq[2];
    mtxq->t[2] = qtrunc(posq[2]);

    //mtx[0][3] = 0;
    //mtx[1][3] = 0;
    //mtx[2][3] = 0;
    //mtx[3][3] = 1;
}

/**
 * Convert float matrix 'src' to fixed point matrix 'dest'.
 * The float matrix may not contain entries larger than 65536 or the console
 * crashes. The fixed point matrix has entries with a 16-bit integer part, so
 * the floating point numbers are multiplied by 2^16 before being cast to a s32
 * integer. If this doesn't fit, the N64 and iQue consoles will throw an
 * exception. On Wii and Wii U Virtual Console the value will simply be clamped
 * and no crashes occur.
 */
void mtxq_to_mtx(Mtx *dest, const ShortMatrix* src) {
    s32 asFixedPoint;
    //register s32 i;
    register s16 *a3 = (s16 *) dest;      // all integer parts stored in first 16 bytes
    register s16 *t0 = (s16 *) dest + 16; // all fraction parts stored in last 16 bytes

    for(int y = 0; y < 3; y++) {
		for(int x = 0; x < 3; x++) {
			asFixedPoint = (s32) src->m[y][x] << 4;
			*a3++ = GET_HIGH_S16_OF_32(asFixedPoint); // integer part
			*t0++ = GET_LOW_S16_OF_32(asFixedPoint); // fraction part
		}
		*a3++ = 0; // integer part
		*t0++ = 0; // fraction part
    }
	for(int x = 0; x < 3; x++) {
		asFixedPoint = src->t[x] << 4;
		*a3++ = GET_HIGH_S16_OF_32(asFixedPoint); // integer part
		*t0++ = GET_LOW_S16_OF_32(asFixedPoint); // fraction part
	}
}

/**
 * Extract a position given an object's transformation matrix and a camera matrix.
 * This is used for determining the world position of the held object: since objMtx
 * inherits the transformation from both the camera and Mario, it calculates this
 * by taking the camera matrix and inverting its transformation by first rotating
 * objMtx back from screen orientation to world orientation, and then subtracting
 * the camera position.
 */
void get_pos_from_transform_mtxq(Vec3q destq, const ShortMatrix* objMtxq, const ShortMatrix* camMtxq) {
	// TODO: optimize with GTE
    q32 camXq = camMtxq->t[0] * camMtxq->m[0][0] + camMtxq->t[1] * camMtxq->m[0][1] + camMtxq->t[2] * camMtxq->m[0][2];
    q32 camYq = camMtxq->t[0] * camMtxq->m[1][0] + camMtxq->t[1] * camMtxq->m[1][1] + camMtxq->t[2] * camMtxq->m[1][2];
    q32 camZq = camMtxq->t[0] * camMtxq->m[2][0] + camMtxq->t[1] * camMtxq->m[2][1] + camMtxq->t[2] * camMtxq->m[2][2];

    destq[0] = objMtxq->t[0] * camMtxq->m[0][0] + objMtxq->t[1] * camMtxq->m[0][1] + objMtxq->t[2] * camMtxq->m[0][2] - camXq;
    destq[1] = objMtxq->t[0] * camMtxq->m[1][0] + objMtxq->t[1] * camMtxq->m[1][1] + objMtxq->t[2] * camMtxq->m[1][2] - camYq;
    destq[2] = objMtxq->t[0] * camMtxq->m[2][0] + objMtxq->t[1] * camMtxq->m[2][1] + objMtxq->t[2] * camMtxq->m[2][2] - camZq;
}

/**
 * Take the vector starting at 'from' pointed at 'to' an retrieve the length
 * of that vector, as well as the yaw and pitch angles.
 * Basically it converts the direction to spherical coordinates.
 */
void vec3f_get_dist_and_angle(Vec3f from, Vec3f to, f32 *dist, s16 *pitch, s16 *yaw) {
    register f32 x = to[0] - from[0];
    register f32 y = to[1] - from[1];
    register f32 z = to[2] - from[2];

    *dist = sqrtf(x * x + y * y + z * z);
    *pitch = atan2s(sqrtf(x * x + z * z), y);
    *yaw = atan2s(z, x);
}

void vec3q_get_dist_and_angle(Vec3q fromq, Vec3q toq, q32 *distq, s16 *pitch, s16 *yaw) {
#ifdef USE_FLOATS
	Vec3f fromf, tof;
	vec3q_to_vec3f(fromf, fromq);
	vec3q_to_vec3f(tof, toq);
	f32 distf = qtof(*distq);
	vec3f_get_dist_and_angle(fromf, tof, &distf, pitch, yaw);
	*distq = q(distf);
#else
    register s32 x = qtrunc(toq[0] - fromq[0]);
    register s32 y = qtrunc(toq[1] - fromq[1]);
    register s32 z = qtrunc(toq[2] - fromq[2]);

	s64 xz = (s64) x * x + (s64) z * z;
    *distq = q(sqrtu(xz + (s64) y * y));
    *pitch = atan2sq(sqrtu((s32) xz), y);
    *yaw = atan2sq(z, x);
#endif
}

/**
 * Construct the 'to' point which is distance 'dist' away from the 'from' position,
 * and has the angles pitch and yaw.
 */
void vec3f_set_dist_and_angle(Vec3f from, Vec3f to, f32 dist, s16 pitch, s16 yaw) {
    to[0] = from[0] + dist * coss(pitch) * sins(yaw);
    to[1] = from[1] + dist * sins(pitch);
    to[2] = from[2] + dist * coss(pitch) * coss(yaw);
}
void vec3q_set_dist_and_angle(Vec3q fromq, Vec3q toq, q32 distq, s16 pitch, s16 yaw) {
#ifdef USE_FLOATS
	Vec3f fromf, tof;
	vec3q_to_vec3f(fromf, fromq);
	vec3q_to_vec3f(tof, toq);
	vec3f_set_dist_and_angle(fromf, tof, qtof(distq), pitch, yaw);
	vec3f_to_vec3q(fromq, fromf);
	vec3f_to_vec3q(toq, tof);
#else
	q32 pitchsinq = sinqs(pitch);
	q32 distpitchcosq = qmul(distq, cosqs(pitch));
	q32 yawsinq = sinqs(yaw);
	q32 yawcosq = cosqs(yaw);
    toq[0] = fromq[0] + qmul(distpitchcosq, yawsinq);
    toq[1] = fromq[1] + qmul(distq, pitchsinq);
    toq[2] = fromq[2] + qmul(distpitchcosq, yawcosq);
#endif
}

/**
 * Return the value 'current' after it tries to approach target, going up at
 * most 'inc' and going down at most 'dec'.
 */
s32 approach_s32(s32 current, s32 target, s32 inc, s32 dec) {
    //! If target is close to the max or min s32, then it's possible to overflow
    // past it without stopping.

    if (current < target) {
        current += inc;
        if (current > target) {
            current = target;
        }
    } else {
        current -= dec;
        if (current < target) {
            current = target;
        }
    }
    return current;
}

/**
 * Return the value 'current' after it tries to approach target, going up at
 * most 'inc' and going down at most 'dec'.
 */
f32 approach_f32(f32 current, f32 target, f32 inc, f32 dec) {
    if (current < target) {
        current += inc;
        if (current > target) {
            current = target;
        }
    } else {
        current -= dec;
        if (current < target) {
            current = target;
        }
    }
    return current;
}

/**
 * Helper function for atan2s. Does a look up of the arctangent of y/x assuming
 * the resulting angle is in range [0, 0x2000] (1/8 of a circle).
 */
static u16 atan2_lookup(f32 y, f32 x) {
    u16 ret;

    if (x == 0) {
        ret = gArctanTable[0];
    } else {
        ret = gArctanTable[(s32)(y / x * 1024 + 0.5f)];
    }
    return ret;
}

static u16 atan2_lookupq(q32 yq, q32 xq) {
#ifdef USE_FLOATS
	return atan2_lookup(qtof(yq), qtof(xq));
#else
    int idx;
    if (xq == 0) {
        idx = 0;
    } else {
        idx = (q64) yq * 1024 / xq;
    }
    return gArctanTable[idx];
#endif
}

/**
 * Compute the angle from (0, 0) to (x, y) as a s16. Given that terrain is in
 * the xz-plane, this is commonly called with (z, x) to get a yaw angle.
 */
s16 atan2s(f32 y, f32 x) {
    u16 ret;

    if (x >= 0) {
        if (y >= 0) {
            if (y >= x) {
                ret = atan2_lookup(x, y);
            } else {
                ret = 0x4000 - atan2_lookup(y, x);
            }
        } else {
            y = -y;
            if (y < x) {
                ret = 0x4000 + atan2_lookup(y, x);
            } else {
                ret = 0x8000 - atan2_lookup(x, y);
            }
        }
    } else {
        x = -x;
        if (y < 0) {
            y = -y;
            if (y >= x) {
                ret = 0x8000 + atan2_lookup(x, y);
            } else {
                ret = 0xC000 - atan2_lookup(y, x);
            }
        } else {
            if (y < x) {
                ret = 0xC000 + atan2_lookup(y, x);
            } else {
                ret = -atan2_lookup(x, y);
            }
        }
    }
    return ret;
}

s16 atan2sq(q32 yq, q32 xq) {
    u16 ret;
    if (xq >= 0) {
        if (yq >= 0) {
            if (yq >= xq) {
                ret = atan2_lookupq(xq, yq);
            } else {
                ret = 0x4000 - atan2_lookupq(yq, xq);
            }
        } else {
            yq = -yq;
            if (yq < xq) {
                ret = 0x4000 + atan2_lookupq(yq, xq);
            } else {
                ret = 0x8000 - atan2_lookupq(xq, yq);
            }
        }
    } else {
        xq = -xq;
        if (yq < 0) {
            yq = -yq;
            if (yq >= xq) {
                ret = 0x8000 + atan2_lookupq(xq, yq);
            } else {
                ret = 0xC000 - atan2_lookupq(yq, xq);
            }
        } else {
            if (yq < xq) {
                ret = 0xC000 + atan2_lookupq(yq, xq);
            } else {
                ret = -atan2_lookupq(xq, yq);
            }
        }
    }
    return ret;
}

/**
 * Compute the atan2 in radians by calling atan2s and converting the result.
 */
f32 atan2f(f32 y, f32 x) {
    return (f32) atan2s(y, x) * M_PI / 0x8000;
}

#define CURVE_BEGIN_1 1
#define CURVE_BEGIN_2 2
#define CURVE_MIDDLE 3
#define CURVE_END_1 4
#define CURVE_END_2 5

/**
 * Set 'result' to a 4-vector with weights corresponding to interpolation
 * value t in [0, 1] and gSplineState. Given the current control point P, these
 * weights are for P[0], P[1], P[2] and P[3] to obtain an interpolated point.
 * The weights naturally sum to 1, and they are also always in range [0, 1] so
 * the interpolated point will never overshoot. The curve is guaranteed to go
 * through the first and last point, but not through intermediate points.
 *
 * gSplineState ensures that the curve is clamped: the first two points
 * and last two points have different weight formulas. These are the weights
 * just before gSplineState transitions:
 * 1: [1, 0, 0, 0]
 * 1->2: [0, 3/12, 7/12, 2/12]
 * 2->3: [0, 1/6, 4/6, 1/6]
 * 3->3: [0, 1/6, 4/6, 1/6] (repeats)
 * 3->4: [0, 1/6, 4/6, 1/6]
 * 4->5: [0, 2/12, 7/12, 3/12]
 * 5: [0, 0, 0, 1]
 *
 * I suspect that the weight formulas will give a 3rd degree B-spline with the
 * common uniform clamped knot vector, e.g. for n points:
 * [0, 0, 0, 0, 1, 2, ... n-1, n, n, n, n]
 * TODO: verify the classification of the spline / figure out how polynomials were computed
 */
void spline_get_weights(Vec4f result, f32 t, UNUSED s32 c) {
    f32 tinv = 1 - t;
    f32 tinv2 = tinv * tinv;
    f32 tinv3 = tinv2 * tinv;
    f32 t2 = t * t;
    f32 t3 = t2 * t;

    switch (gSplineState) {
        case CURVE_BEGIN_1:
            result[0] = tinv3;
            result[1] = t3 * 1.75f - t2 * 4.5f + t * 3.0f;
            result[2] = -t3 * (11 / 12.0f) + t2 * 1.5f;
            result[3] = t3 * (1 / 6.0f);
            break;
        case CURVE_BEGIN_2:
            result[0] = tinv3 * 0.25f;
            result[1] = t3 * (7 / 12.0f) - t2 * 1.25f + t * 0.25f + (7 / 12.0f);
            result[2] = -t3 * 0.5f + t2 * 0.5f + t * 0.5f + (1 / 6.0f);
            result[3] = t3 * (1 / 6.0f);
            break;
        case CURVE_MIDDLE:
            result[0] = tinv3 * (1 / 6.0f);
            result[1] = t3 * 0.5f - t2 + (4 / 6.0f);
            result[2] = -t3 * 0.5f + t2 * 0.5f + t * 0.5f + (1 / 6.0f);
            result[3] = t3 * (1 / 6.0f);
            break;
        case CURVE_END_1:
            result[0] = tinv3 * (1 / 6.0f);
            result[1] = -tinv3 * 0.5f + tinv2 * 0.5f + tinv * 0.5f + (1 / 6.0f);
            result[2] = tinv3 * (7 / 12.0f) - tinv2 * 1.25f + tinv * 0.25f + (7 / 12.0f);
            result[3] = t3 * 0.25f;
            break;
        case CURVE_END_2:
            result[0] = tinv3 * (1 / 6.0f);
            result[1] = -tinv3 * (11 / 12.0f) + tinv2 * 1.5f;
            result[2] = tinv3 * 1.75f - tinv2 * 4.5f + tinv * 3.0f;
            result[3] = t3;
            break;
    }
}

/**
 * Initialize a spline animation.
 * 'keyFrames' should be an array of (s, x, y, z) vectors
 *  s: the speed of the keyframe in 1000/frames, e.g. s=100 means the keyframe lasts 10 frames
 *  (x, y, z): point in 3D space on the curve
 * The array should end with three entries with s=0 (infinite keyframe duration).
 * That's because the spline has a 3rd degree polynomial, so it looks 3 points ahead.
 */
void anim_spline_init(Vec4s *keyFrames) {
    gSplineKeyframe = keyFrames;
    gSplineKeyframeFraction = 0;
    gSplineState = 1;
}

/**
 * Poll the next point from a spline animation.
 * anim_spline_init should be called before polling for vectors.
 * Returns TRUE when the last point is reached, FALSE otherwise.
 */
s32 anim_spline_poll(Vec3f result) {
    Vec4f weights;
    s32 i;
    s32 hasEnded = FALSE;

    for(int i = 0; i < 3; i++) result[i] = 0;
    spline_get_weights(weights, gSplineKeyframeFraction, gSplineState);
    for (i = 0; i < 4; i++) {
        result[0] += weights[i] * gSplineKeyframe[i][1];
        result[1] += weights[i] * gSplineKeyframe[i][2];
        result[2] += weights[i] * gSplineKeyframe[i][3];
    }

    if ((gSplineKeyframeFraction += gSplineKeyframe[0][0] / 1000.0f) >= 1) {
        gSplineKeyframe++;
        gSplineKeyframeFraction--;
        switch (gSplineState) {
            case CURVE_END_2:
                hasEnded = TRUE;
                break;
            case CURVE_MIDDLE:
                if (gSplineKeyframe[2][0] == 0) {
                    gSplineState = CURVE_END_1;
                }
                break;
            default:
                gSplineState++;
                break;
        }
    }

    return hasEnded;
}
