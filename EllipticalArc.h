// Copyright (c) 2003-2004, Luc Maisonobe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that
// the following conditions are met:
//
//    Redistributions of source code must retain the
//    above copyright notice, this list of conditions and
//    the following disclaimer.
//    Redistributions in binary form must reproduce the
//    above copyright notice, this list of conditions and
//    the following disclaimer in the documentation
//    and/or other materials provided with the
//    distribution.
//    Neither the names of spaceroots.org, spaceroots.com
//    nor the names of their contributors may be used to
//    endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
// CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Ported to Objective-C for CoreGraphics on 21. August 2012 by Niklas Hauser

#import <Foundation/Foundation.h>

/** This class represents an elliptical arc on a 2D plane.
 *
 * <p>This class can handle parts of ellipse in addition to full
 * ellipses and it can handle ellipses which are not aligned with the
 * x and y reference axes of the plane.</p>
 *
 * <p>Another improvement is that this class can handle degenerated
 * cases like for example very flat ellipses (semi-minor axis much
 * smaller than semi-major axis) and drawing of very small parts of
 * such ellipses at very high magnification scales. This imply
 * monitoring the drawing approximation error for extremely small
 * values. Such cases occur for example while drawing orbits of comets
 * near the perihelion.</p>
 *
 * <p>When the arc does not cover the complete ellipse, the lines
 * joining the center of the ellipse to the endpoints can optionally
 * be included or not in the outline, hence allowing to use it for
 * pie-charts rendering. If these lines are not included, the curve is
 * not naturally closed.</p>
 *
 * @author L. Maisonobe
 * @author Niklas Hauser
 */
@interface EllipticalArc : NSObject
{
	/** Abscissa of the center of the ellipse. */
	double _cx;
	/** Ordinate of the center of the ellipse. */
	double _cy;
	/** Semi-major axis. */
	double _a;
	/** Semi-minor axis. */
	double _b;
	/** Orientation of the major axis with respect to the x axis. */
	double _theta;
	double _cosTheta;
	double _sinTheta;
	/** Start angle of the arc. */
	double _eta1;
	/** End angle of the arc. */
	double _eta2;
	/** Abscissa of the start point. */
	double _x1;
	/** Ordinate of the start point. */
	double _y1;
	/** Abscissa of the end point. */
	double _x2;
	/** Ordinate of the end point. */
	double _y2;
	/** Abscissa of the first focus. */
	double _xF1;
	/** Ordinate of the first focus. */
	double _yF1;
	/** Abscissa of the second focus. */
	double _xF2;
	/** Ordinate of the second focus. */
	double _yF2;
	/** Abscissa of the leftmost point of the arc. */
	double _xLeft;
	/** Ordinate of the highest point of the arc. */
	double _yUp;
	/** Horizontal width of the arc. */
	double _width;
	/** Vertical height of the arc. */
	double _height;
	/** Indicator for center to endpoints line inclusion. */
	BOOL _isPieSlice;
	/** Maximal degree for Bézier curve approximation. */
	int _maxDegree;
	/** Default flatness for Bézier curve approximation. */
	double _defaultFlatness;
	
	double _f;
	double _e2;
	double _g;
	double _g2;
}

/** Simple constructor.
 * Build an elliptical arc composed of the full unit circle centered
 * on origin
 */
-(id)init;

/** Build an elliptical arc from its canonical geometrical elements.
 * @param center center of the ellipse
 * @param a semi-major axis
 * @param b semi-minor axis
 * @param theta orientation of the major axis with respect to the x axis
 * @param lambda1 start angle of the arc
 * @param lambda2 end angle of the arc
 * @param pieSlice if true, the lines between the center of the ellipse
 * and the endpoints are part of the shape (it is pie slice like)
 */
-(id)initWithCentre:(CGPoint)centre
	  semiMajorAxis:(double)a
	  semiMinorAxis:(double)b
		orientation:(double)theta
		 startAngle:(double)lambda1
		   endAngle:(double)lambda2
		   pieSlice:(BOOL)pieSlice;

/** Build a full ellipse from its canonical geometrical elements.
 * @param center center of the ellipse
 * @param a semi-major axis
 * @param b semi-minor axis
 * @param theta orientation of the major axis with respect to the x axis
 */
-(id)initWithCentre:(CGPoint)centre
	  semiMajorAxis:(double)a
	  semiMinorAxis:(double)b
		orientation:(double)theta;

/** Get the elliptical arc point for a given angular parameter.
 * @param lambda angular parameter for which point is desired
 * @return the desired elliptical arc point location
 */
-(CGPoint)pointAtAngle:(double)lambda;

/** Tests if the specified coordinates are inside the boundary of the Shape.
 * @param point test point
 * @return true if the specified coordinates are inside the Shape
 * boundary; false otherwise
 */
-(BOOL)containsPoint:(CGPoint)point;

/** Tests if the interior of the Shape entirely contains the
 * specified rectangular area.
 * @param rect rectangle area
 * @return true if the interior of the Shape entirely contains the
 * specified rectangular area; false otherwise
 */
-(BOOL)containsRectangle:(CGRect)rect;

/** Returns a bounding box rectangle that completely encloses the Shape.
 */
-(CGRect)bounds;

/** Tests if the interior of the Shape intersects the interior of a
 * specified rectangle.
 */
-(BOOL)intersectsRectangle:(CGRect)rect;

/** 
 * @param m affine transformation, set to NULL if not desired
 * @return a transformed CGMutablePathRef to use with CoreGrapics.
 */
-(CGMutablePathRef)pathWithTranformation:(CGAffineTransform *)m;

/**
 * @param m affine transformation, set to NULL if not desired
 * @param flatness desired maximum error
 * @return a transformed CGMutablePathRef to use with CoreGrapics.
 */
-(CGMutablePathRef)pathWithTranformation:(CGAffineTransform *)m flatness:(double)flatness;

@end
