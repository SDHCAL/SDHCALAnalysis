  /// \file CutProcessor.cc
/*
 *
 * CutProcessor.cc source template automatically generated by a class generator
 * Creation date : lun. janv. 5 2015
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


#include "Processors/CutProcessor.h"

#include "marlin/Global.h"
#include "marlin/VerbosityLevels.h"

// lcio includes
#include "UTIL/CellIDDecoder.h"
#include "EVENT/CalorimeterHit.h"

// std
#include <limits>
#include <cstdio>

// sdhcal
#include "Common/Linear3DFit.h"

CutProcessor aCutProcessor;

CutProcessor::CutProcessor() :
	marlin::Processor("CutProcessor")
{
	_description = "CutProcessor to apply some basic cuts on sdhcal data. Cut electrons and muons, keep only pions";

	registerProcessorParameter("NHitCut",
				 "cut on total number of hits",
				 m_nHitCut,
				 static_cast<int> (100) );

	registerProcessorParameter("LayerFractionCut",
				 "cut on nTouchedLayer/NLayer",
				 m_layerFractionCut,
				 static_cast<double> (0.2) );

	registerProcessorParameter("NHitOverNLayerCut",
				 "cut on NHit/nbOfTouchedLayers",
				 m_nHitOverNLayerCut,
				 static_cast<double> (3) );

	registerProcessorParameter("RadiusOverCog2Cut",
				 "cut on radius/cogz",
				 m_radiusOverCog2Cut,
				 static_cast<double> (0.4) );

	registerProcessorParameter("ShowerStartingLayerCut",
				 "cut on shower starting layer",
				 m_showerStartingLayerCut,
				 static_cast<int> (0) );

	registerProcessorParameter("NTouchedLayersCut",
				 "cut on the number of touched layers",
				 m_nTouchedLayersCut,
				 static_cast<int> (0) );

	registerProcessorParameter("FractalTimesCentralCellsCut",
				 "cut on (fractalDim * hitsInCentralCells) / (ln(NHit) * NHit)",
				 m_fractalTimesCentralCellsCut,
				 static_cast<double> (0.0) );

	registerProcessorParameter("NHolesCut",
				 "cut on the number of holes after the shower starting point",
				 m_nHolesCut,
				 static_cast<int> (0) );

	registerProcessorParameter("NHitEdgePercentCut",
				 "cut on the percentage of hits in the edges",
				 m_nHitEdgePercentCut,
				 static_cast<double> (0.5) );

	registerProcessorParameter("LargeRMSCut",
				 "cut on rms of the 5 first layers",
				 m_largeRMSCut,
				 static_cast<double> (4.0) );

	registerProcessorParameter("BarycenterPositionCut",
				 "cut on barcenter position in x-y (central region)",
				 m_barycenterPositionCut,
				 static_cast<double> (4.0) );

	registerProcessorParameter("CosThetaCut",
				 "cut on the cos theta of the shower fit",
				 m_cosThetaCut,
				 static_cast<double> (0.95) );

	registerProcessorParameter("NeutralFirstLayerCut",
				 "cut on neutral particles (first layers are empties)",
				 m_neutralFirstLayerCut,
				 static_cast<int> (5) );

	registerProcessorParameter("DecoderString" ,
				 "decoder string for cell ID decoder" ,
				 m_decoderString,
				 std::string("M:3,S-1:3,I:9,J:9,K-1:6"));

	registerProcessorParameter("SDHCALCollectionName" ,
				 "collection name for SDHCAL hits" ,
				 m_sdhcalCollectionName,
				 std::string("HCALBarrel"));

	std::vector<std::string> ijkVec;
	ijkVec.push_back("I");
	ijkVec.push_back("J");
	ijkVec.push_back("K-1");

	registerProcessorParameter("IJKEncoding",
				 "I J K hit encoding",
				 m_ijkEncoding,
				 ijkVec);

}

CutProcessor::~CutProcessor()
{
	/* nop */
}


int CutProcessor::ijkToKey(const int i , const int j , const int k)
{
	return 100*100*k+100*j+i;
}

std::vector<int> CutProcessor::keyToIJK( const int &key )
{

	std::vector<int> vec;
	vec.push_back( key%100 );
	vec.push_back( key/100%100 );
	vec.push_back( key/10000 );
	return vec;
}


bool CutProcessor::sortByLayer( EVENT::CalorimeterHit *caloHit1 , EVENT::CalorimeterHit *caloHit2 )
{
	return ( caloHit1->getPosition()[2] < caloHit2->getPosition()[2] );
}


void CutProcessor::init()
{
	m_nbOfLayers = 48;
	m_nCells0 = 96;
	m_nCells1 = 96;
	m_nMultiParticleEvents = 0;
	m_nProcessedEvents = 0;

	printParameters();
}


void CutProcessor::processRunHeader( LCRunHeader* run )
{
	/* nop */
}


void CutProcessor::processEvent( EVENT::LCEvent * evt )
{
	streamlog_out(DEBUG) << "CutProcessor - processing event no " << evt->getEventNumber() << std::endl;

	m_nProcessedEvents++;
	reset();

	UTIL::CellIDDecoder<CalorimeterHit>::setDefaultEncoding(m_decoderString);

	EVENT::LCCollection *pCollection = 0;

	try
	{
		pCollection = evt->getCollection( m_sdhcalCollectionName );
	}
	catch( DataNotAvailableException &e )
	{
		streamlog_out(ERROR) << "LCIO exception thrown : " << e.what() << std::endl;
		throw marlin::StopProcessingException(this);
	}

	// first cut on nHit
	if(pCollection->getNumberOfElements() < m_nHitCut)
	{
		streamlog_out(DEBUG) << "Skipping event - NHit cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	int nHit1 = 0;
	int nHit2 = 0;
	int nHit3 = 0;

	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double weight = 1.0;
	double sumweight = 0.0;

	std::vector<ThreeVector> positions;

	UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder(pCollection);

	for(unsigned int i=0 ; i<pCollection->getNumberOfElements() ; i++)
	{
		if( pCollection->getTypeName() != LCIO::CALORIMETERHIT )
			throw marlin::StopProcessingException(this);

		CalorimeterHit *caloHit = static_cast<CalorimeterHit*> (pCollection->getElementAt(i));

		float fThr = caloHit->getEnergy();
		int   I = cellIdDecoder(caloHit)[m_ijkEncoding.at(0).c_str()];
		int   J = cellIdDecoder(caloHit)[m_ijkEncoding.at(1).c_str()];
		int   K = cellIdDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		if(K > m_nbOfLayers - 1)
			continue;

		if( I < 5 || I > 95 || J < 5 || J > 95 )
			m_nHitsInEdge ++;

		// the touched layers
		m_nTouchedLayers.at( K ) ++;

		// number of hits for each threshold + weight for cog computation
		if( (fThr-1.f) < std::numeric_limits<float>::epsilon() )
		{
			weight = 1.0;
			nHit1++;
		}
		else if( (fThr-2.0)< std::numeric_limits<float>::epsilon() )
		{
			weight = 5.0;
			nHit2++;
		}
		else if( (fThr-3.0)< std::numeric_limits<float>::epsilon() )
		{
			weight = 10.0;
			nHit3++;
		}


		x += weight*I;
		y += weight*J;
		z += weight*K;
		sumweight += weight;

		m_caloHitCollection.push_back(caloHit);
		positions.push_back(ThreeVector(caloHit->getPosition()[0], caloHit->getPosition()[1], caloHit->getPosition()[2]));
	}

	std::sort(m_caloHitCollection.begin(), m_caloHitCollection.end(), CutProcessor::sortByLayer);

	int nTouchedLayers = 0;
	// get the number of touched layers in the event
	for(unsigned int l=0 ; l<m_nTouchedLayers.size() ; l++)
		if(m_nTouchedLayers.at(l) != 0)
			nTouchedLayers++;

	m_nHit.push_back(m_caloHitCollection.size());
	m_nHit.push_back(nHit1);
	m_nHit.push_back(nHit2);
	m_nHit.push_back(nHit3);

	m_cog.push_back( x/sumweight );
	m_cog.push_back( y/sumweight );
	m_cog.push_back( z/sumweight );

//	m_fractalDimension = this->getFractalDimension();

	// calculate the shower starting layer
	for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ )
	{
		CalorimeterHit *caloHit = m_caloHitCollection.at(h);

		float fThr = caloHit->getEnergy();
		int I = cellIdDecoder(caloHit)[m_ijkEncoding.at(0).c_str()];
		int J = cellIdDecoder(caloHit)[m_ijkEncoding.at(1).c_str()];
		int K = cellIdDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		if(fabs( I - m_cog.at(0) ) > 5 || fabs( J - m_cog.at(1) ) > 5)
			continue;

		int count = 0;
		std::vector<int> count2(3,0);

		for( unsigned int h2=0 ; h2<m_caloHitCollection.size() ; h2++ )
		{
			if( h == h2 )
				continue;

			CalorimeterHit *caloHit2 = m_caloHitCollection.at(h2);

			int I2 = cellIdDecoder(caloHit2)[m_ijkEncoding.at(0).c_str()];
			int J2 = cellIdDecoder(caloHit2)[m_ijkEncoding.at(1).c_str()];
			int K2 = cellIdDecoder(caloHit2)[m_ijkEncoding.at(2).c_str()];

			if( fabs( I2 - m_cog.at(0) ) < 5 && fabs( J2 - m_cog.at(1) ) < 5 )
				count++;
			else
				continue;

			if( K2 == K+1 ) count2.at(0) ++;
			if( K2 == K+2 ) count2.at(1) ++;
			if( K2 == K+3 ) count2.at(2) ++;
		}

		if( count <= 4 )
			continue;

		if( count2.at(0) >= 4 && count2.at(1) >= 4 && count2.at(2) >= 4 )
		{
			m_showerStartingLayer = K;
			break;
		}
	}

	sumweight = 0.0;
	int count = 0;

	for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ )
	{
		CalorimeterHit *caloHit = m_caloHitCollection.at(h);

		int I = cellIdDecoder(caloHit)[m_ijkEncoding.at(0).c_str()];
		int J = cellIdDecoder(caloHit)[m_ijkEncoding.at(1).c_str()];
		int K = cellIdDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		if( K >= m_showerStartingLayer )
		{
			m_radius += (m_cog.at(0) - I)*(m_cog.at(0) - I) + (m_cog.at(1) - J)*(m_cog.at(1) - J);
			sumweight ++;
		}

		if( ( I - m_cog.at(0) ) < 3 && ( J - m_cog.at(1) ) < 3 )
			count ++;
	}

	m_nHitsInCentralCells = count;

	if( sumweight != 0 )
		sqrt( m_radius /= sumweight );


	if(m_showerStartingLayer != -1)
	{
		for( int Kiter = m_showerStartingLayer+1 ; Kiter < m_showerStartingLayer+8  ; Kiter++ ) {

			int count = 0;
			for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ ) {

				CalorimeterHit *caloHit = m_caloHitCollection.at(h);

				float fThr = caloHit->getEnergy();
				int I = cellIdDecoder(caloHit)[m_ijkEncoding.at(0).c_str()];
				int J = cellIdDecoder(caloHit)[m_ijkEncoding.at(1).c_str()];
				int K = cellIdDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

				if( K == Kiter && ( I - m_cog.at(0) ) < 10 && ( J - m_cog.at(1) ) < 10 )
				{
					count++;
					break;
				}
			}
			if( count == 0 )
				m_nHolesAfterStartingPoint++;
		}
	}

	bool isNeutral = true;

	for(unsigned int k=0 ; k<m_neutralFirstLayerCut ; k++)
	{
		if(m_nTouchedLayers.at(k) != 0)
		{
			isNeutral = false;
			break;
		}
	}

	if(isNeutral)
	{
		streamlog_out(DEBUG)  << "Skipping event - Neutral particle cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	if((double)m_nHit.at(0)/(double)nTouchedLayers < m_nHitOverNLayerCut)
	{
		streamlog_out(DEBUG)  << "Skipping event - NHit/NLayers cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	if((double)nTouchedLayers/(double)m_nbOfLayers < m_layerFractionCut)
	{
		streamlog_out(DEBUG)  << "Skipping event - nTouchedLayers/nLayers cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	if(m_radius / m_cog.at(2) < m_radiusOverCog2Cut)
	{
		streamlog_out(DEBUG)  << "Skipping event - radius/cog[2] cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	if(m_showerStartingLayer < m_showerStartingLayerCut && nTouchedLayers < m_nTouchedLayersCut)
	{
		streamlog_out(DEBUG)  << "Skipping event - startPoint && NLayer cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	if((m_fractalDimension * m_nHitsInCentralCells) / (std::log(m_nHit.at(0)) * m_nHit.at(0)) < m_fractalTimesCentralCellsCut)
	{
		streamlog_out(DEBUG)  << "Skipping event - fractal dim cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	if(m_nHolesAfterStartingPoint > m_nHolesCut)
	{
		streamlog_out(DEBUG)  << "Skipping event - NHoles cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	if(m_nHitsInEdge / m_nHit.at(0) > m_nHitEdgePercentCut)
	{
		streamlog_out(DEBUG)  << "Skipping event - NEdge/NHit cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	// cog too fare from center, then cut
	if(m_cog.at(0) < m_barycenterPositionCut || m_nCells0 - m_cog.at(0) < m_barycenterPositionCut
	|| m_cog.at(1) < m_barycenterPositionCut || m_nCells1 - m_cog.at(1) < m_barycenterPositionCut)
	{
		streamlog_out(DEBUG)  << "Skipping event - Cog border cut" << std::endl;
		throw marlin::SkipEventException(this);
	}

	Linear3DFit *pFitter = new Linear3DFit(positions);
	pFitter->fit();

	ThreeVector px(-1, 0, pFitter->getFitParameters()[1]);
	ThreeVector py(0, -1, pFitter->getFitParameters()[3]);
	float cosTheta = px.cross(py).cosTheta();

	delete pFitter;

	if(cosTheta < m_cosThetaCut)
	{
		streamlog_out(DEBUG)  << "Skipping event - cos theta shower cut (" << cosTheta << " < " << m_cosThetaCut << ")" << std::endl;
		throw marlin::SkipEventException(this);
	}

	bool isSingleParticle = this->isSingleParticle();
	std::string single = isSingleParticle ? "true" : "false";

	if(!isSingleParticle)
	{
		m_nMultiParticleEvents++;
		m_evtIdMultiEvent.push_back(evt->getEventNumber());

		streamlog_out(DEBUG)  << "Skipping event - multi event cut" << std::endl;
		throw marlin::SkipEventException(this);
	}
}

//---------------------------------------------------------------------------------------------

void CutProcessor::check( LCEvent *evt )
{
	/* nop */
}

//---------------------------------------------------------------------------------------------

void CutProcessor::end()
{
	reset();
	streamlog_out(DEBUG)  << "Number of multi-particle events : " << m_nMultiParticleEvents << " , fraction = " << (double(m_nMultiParticleEvents)/double(m_nProcessedEvents))*100 << std::endl;

	streamlog_out(DEBUG) << "Event id with event : " << std::endl;
	for(unsigned int i=0 ; i<m_evtIdMultiEvent.size() ; i++)
	{
		streamlog_out(DEBUG) << m_evtIdMultiEvent.at(i) << std::endl;
	}
}

//---------------------------------------------------------------------------------------------

bool CutProcessor::isSingleParticle()
{
	std::vector<double> layerCogX(m_nbOfLayers, 0);
	std::vector<double> layerCogY(m_nbOfLayers, 0);
	std::vector<double> layerRMSX(m_nbOfLayers, 0);
	std::vector<double> layerRMSY(m_nbOfLayers, 0);
	std::vector<double> weightSum(m_nbOfLayers, 0);

	UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(m_decoderString);

	int firstLayersCut = 5;
	double largeRMSCut = 5.0;
	int largeRMSCounterCut = 4;
	std::vector<EVENT::CalorimeterHit*> firstCaloHitCollection;

	for(unsigned int i=0 ; i<m_caloHitCollection.size() ; i++)
	{
		CalorimeterHit *caloHit = m_caloHitCollection.at(i);
		int K = idDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		if(K < firstLayersCut)
			firstCaloHitCollection.push_back(caloHit);
	}

	// compute the cog for each layer
	for(unsigned int i=0 ; i<firstCaloHitCollection.size() ; i++)
	{
		CalorimeterHit *caloHit = firstCaloHitCollection.at(i);

		float fThr = caloHit->getEnergy();

		int I = idDecoder(caloHit)[m_ijkEncoding.at(0).c_str()];
		int J = idDecoder(caloHit)[m_ijkEncoding.at(1).c_str()];
		int K = idDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		double weight = 0.;
		if( (fThr-1.f) < std::numeric_limits<float>::epsilon() )
			weight = 1.0;
		else if( (fThr-2.0)< std::numeric_limits<float>::epsilon() )
			weight = 5.0;
		else if( (fThr-3.0)< std::numeric_limits<float>::epsilon() )
			weight = 10.0;

		layerCogX.at(K) += I*weight;
		layerCogY.at(K) += J*weight;
		weightSum.at(K) += weight;
	}

	// normalize the cog
	for(unsigned int i=0 ; i<layerCogX.size() ; i++)
	{
		if(weightSum.at(i) != 0)
		{
			layerCogX.at(i) /= double(weightSum.at(i));
			layerCogY.at(i) /= double(weightSum.at(i));
		}
	}


	// compute the rms for each layer
	for(unsigned int i=0 ; i<firstCaloHitCollection.size() ; i++)
	{
		CalorimeterHit *caloHit = firstCaloHitCollection.at(i);

		float fThr = caloHit->getEnergy();

		int I = idDecoder(caloHit)[m_ijkEncoding.at(0).c_str()];
		int J = idDecoder(caloHit)[m_ijkEncoding.at(1).c_str()];
		int K = idDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		double weight = 0.;
		if( (fThr-1.f) < std::numeric_limits<float>::epsilon() )
			weight = 1.0;
		else if( (fThr-2.0)< std::numeric_limits<float>::epsilon() )
			weight = 5.0;
		else if( (fThr-3.0)< std::numeric_limits<float>::epsilon() )
			weight = 10.0;

		layerRMSX.at(K) += (I - layerCogX.at(K))*(I - layerCogX.at(K))*weight;
		layerRMSY.at(K) += (J - layerCogY.at(K))*(J - layerCogY.at(K))*weight;
	}

	int largeRMSCounter = 0;

	// normalize the rms
	for(unsigned int i=0 ; i<layerCogX.size() ; i++)
	{
		if(weightSum.at(i) != 0)
		{
			layerRMSX.at(i) = std::sqrt(layerRMSX.at(i)/double(weightSum.at(i)));
			layerRMSY.at(i) = std::sqrt(layerRMSY.at(i)/double(weightSum.at(i)));

//			streamlog_out(DEBUG) << "RMS in layer " << i << " is (" << layerRMSX.at(i) << " , " << layerRMSY.at(i) << ")" << std::endl;

			if(layerRMSX.at(i) > m_largeRMSCut || layerRMSY.at(i) > m_largeRMSCut)
				largeRMSCounter++;
		}
	}

//	streamlog_out(DEBUG) << "largeRMSCounter is " << largeRMSCounter << std::endl;

	return largeRMSCounter < largeRMSCounterCut ? true : false;
}

//---------------------------------------------------------------------------------------------

double CutProcessor::getFractalDimension()
{

	int vec[] = {2,3,4,6,8,12,16};
	float f3D = 0;

	for(int i=0; i<7; i++) {

		int nCube = nHitsInCube(vec[i]);
		if(nCube >= m_nHit.at(0)) return 0;
		f3D += std::log(float(m_nHit.at(0))/nCube)/std::log(vec[i]);
	}
	return f3D/7.0;
}

//---------------------------------------------------------------------------------------------

int CutProcessor::nHitsInCube( int cubeSize )
{
	std::vector<int> keys;
	int ncube = 0;
	UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(m_decoderString);

	for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ ) {

		CalorimeterHit *caloHit = m_caloHitCollection.at(h);

		int newI = idDecoder(caloHit)[m_ijkEncoding.at(0).c_str()]/(cubeSize+1);
		int newJ = idDecoder(caloHit)[m_ijkEncoding.at(1).c_str()]/(cubeSize+1);
		int newK = idDecoder(caloHit)[m_ijkEncoding.at(2).c_str()]/(cubeSize+1);

		int key = ijkToKey( newI , newJ , newK );

		if( std::find( keys.begin() , keys.end() , key ) != keys.end() )
			continue;

		ncube++;
		keys.push_back( key );
	}
	return ncube;
}

//---------------------------------------------------------------------------------------------

void CutProcessor::reset()
{
	m_nTouchedLayers = std::vector<int>(m_nbOfLayers, 0);
	m_radius = 0.0;
	m_showerStartingLayer = -1;
	m_nHitsInEdge = 0;
	m_nHolesAfterStartingPoint = 0;
	m_fractalDimension = 0.0;
	m_nHitsInCentralCells = 0;

	m_nHit.clear();
	m_cog.clear();
	m_caloHitCollection.clear();
}


