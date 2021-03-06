  /// \file Linear3DFit.cc
/*
 *
 * Linear3DFit.cc source template automatically generated by a class generator
 * Creation date : mar. janv. 6 2015
 *
 * This file is part of SDHCAL libraries.
 * 
 * SDHCAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * SDHCAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SDHCAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author rete
 * @copyright CNRS , IPNL
 */


#include "Common/Linear3DFit.h"

Linear3DFit::Linear3DFit(const std::vector<ThreeVector> &pos)
{
	chi2 = params[0] = params[1] = params[2] = params[3] = 0;
	paramsError[0] = paramsError[1] = paramsError[2] = paramsError[3] = 0;
	positions = pos;
}


Linear3DFit::~Linear3DFit()
{

}

void Linear3DFit::computeChi2()
{
	chi2 = 0;

	for( unsigned int i=0 ; i<positions.size() ; i++ )
	{
		ThreeVector temp(this->getNormaleProjection( positions.at(i) ));
		float d = this->getNormaleProjection( positions.at(i) ).mag();
		float err = 1./sqrt(12);
		chi2 += (d/err)*(d/err);
	}

	chi2 = chi2/(positions.size()-1);
}


// The position of the projection of the point on the real line
ThreeVector Linear3DFit::vectorFromRealLine(const ThreeVector& vec) const
{
	ThreeVector x0( params[0], params[2] , 0. );
	ThreeVector x1( params[0] + params[1], params[2] + params[3] , 1. );
	ThreeVector u = (x1-x0);
	u.setMag( (vec-x1).mag() * cos( u.angle( vec-x1 ) ) );

	return (x1 + u);
}

// The vector that join the given point (vec) and his normal projection on the line
ThreeVector Linear3DFit::getNormaleProjection(const ThreeVector& vec) const
{
	return ( vec - this->vectorFromRealLine(vec) );
}

void Linear3DFit::fit()
{
	if(positions.empty())
		return;

	params[0] = params[1] = params[2] = params[3] = 0;
	paramsError[0] = paramsError[1] = paramsError[2] = paramsError[3] = 0;

	float xsum = 0.0;
	float ysum = 0.0;
	float zsum = 0.0;
	float zzsum = 0.0;
	float xzsum = 0.0;
	float yzsum = 0.0;

	for ( unsigned int i=0 ; i<positions.size() ; i++ )
	{
		//for equation 1
		zsum = zsum + positions.at(i).z();
		xsum = xsum + positions.at(i).x();
		zzsum = zzsum + (positions.at(i).z()*positions.at(i).z());
		xzsum = xzsum + positions.at(i).x()*positions.at(i).z();

		//for equation 2
		ysum = ysum + positions.at(i).y();
		yzsum = yzsum + positions.at(i).y()*positions.at(i).z();
	}

	float A1 = zsum;
	float B1 = positions.size();
	float C1 = xsum;
	float D1 = zzsum;
	float E1 = xzsum;
	float C2 = ysum;
	float E2 = yzsum;

	params[0] = (D1*C1-E1*A1)/(B1*D1-A1*A1);
	params[1] = (E1*B1-C1*A1)/(B1*D1-A1*A1);
	params[2] = (D1*C2-E2*A1)/(B1*D1-A1*A1);
	params[3] = (E2*B1-C2*A1)/(B1*D1-A1*A1);

	paramsError[0] = sqrt(D1/(B1*D1-A1*A1));
	paramsError[1] = sqrt(B1/(B1*D1-A1*A1));
	paramsError[2] = sqrt(D1/(B1*D1-A1*A1));
	paramsError[3] = sqrt(B1/(B1*D1-A1*A1));

	computeChi2();
}

void Linear3DFit::fit(const std::vector<ThreeVector> &pos)
{
	positions = pos;
	fit();
}

float Linear3DFit::getChi2() const
{
	return chi2;
}

const float *Linear3DFit::getFitParameters() const
{
	return params;
}

const float *Linear3DFit::getFitParError() const
{
	return paramsError;
}
