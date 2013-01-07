#include "GeometricAnalysis.h"

GeometricAnalysis::GeometricAnalysis ()
	:  SD1(0), SD2(0), TINN(0), M(0), N(0), HRVTriangularIndex(0), Y(0), X(0), HistogramBinLength(7.8125), SamplingInterval(1/360)
{}

GeometricAnalysis::~GeometricAnalysis (void)
{}

void GeometricAnalysis::runModule(const ECGInfo &info, const ECGRs &ecgRs, ECGHRV2 &ecgHRV2)
{
	this->rpeaks = ecgRs;
	this->hrv2 = ecgHRV2;
	this->SamplingInterval = 1/info.channel_one.frequecy; 
	PrepareRRSignal();
	MakeHistogramAndGeometricalParams();
	MakePoincareAndSDParams();
	SetHRV2Params();
}

void GeometricAnalysis::PrepareRRSignal()
{
	OtherSignal RR_intervals;

	auto rpeaks_size = rpeaks.GetRs()->signal->size;
	gsl_vector_int_view RR1 = gsl_vector_int_subvector(rpeaks.GetRs()->signal,1,rpeaks_size-1);
	gsl_vector_int_view RR2 = gsl_vector_int_subvector(rpeaks.GetRs()->signal,0,rpeaks_size-1);
	gsl_vector_int_sub(&RR1.vector,&RR2.vector);

	for(int i = 0; i < rpeaks_size-1; i++) 
	{
		auto Value = gsl_vector_int_get (&RR1.vector, i);
		gsl_vector_set(RR_intervals->signal, i, Value);
	}

	gsl_vector_scale(RR_intervals->signal,SamplingInterval);

	this->RR_intervals=RR_intervals;
}

void GeometricAnalysis::MakeHistogramAndGeometricalParams()
{

	IntSignal Histogram;

	auto RRmax = gsl_vector_max(RR_intervals->signal);
	auto RRmin = gsl_vector_min(RR_intervals->signal);
	auto length = RRmax-RRmin;
	auto index_hist = length/HistogramBinLength;

	auto RR_intervals_size = RR_intervals->signal->size;

	for(int i=0; i < index_hist; i++)
	{
		gsl_vector_int_set(Histogram->signal, i,0);

		for(int j=0; j < RR_intervals_size; j++)
		{
			auto value = gsl_vector_get(RR_intervals->signal, j);

			if(value > (RRmin + (i * HistogramBinLength ) ) )
			{
				if(value < (RRmin + ( (i+1) * HistogramBinLength ) ) )
				{
					auto oldvalue = gsl_vector_int_get(Histogram->signal, i); 
					gsl_vector_int_set(Histogram->signal, i,oldvalue+1);
				}
			}
		}
	}

	auto Y = gsl_vector_int_max(Histogram->signal);
	auto X = gsl_vector_int_max_index(Histogram->signal);

	double minimum = 0;
	double globalminimum = 10000000000000;
	int N = 0;
	int M = 0;
	double x[3] = { N, X, M };
	double y[3] = { 0, Y, 0 };

	gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,3);
	gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

	auto Histogram_size = Histogram->signal->size;

	for(int index_N=0; index_N < X-1; index_N++)
	{
		for(int index_M=X+1; index_M < Histogram_size; index_M++)
		{
			x[0] = index_N; x[2] = index_M;

			gsl_interp_init(interpolation, x, y, 3);

			for(int i=0; i < index_N; i++)
			{
				auto HistogramValue = gsl_vector_int_get(Histogram->signal, i);
				minimum = minimum + (HistogramValue*HistogramValue);
			}

			for(int i=index_N+1; i < index_M-1; i++)
			{
				auto LinearValue = gsl_interp_eval(interpolation, x, y, i, accelerator);
				auto HistogramValue = gsl_vector_int_get(Histogram->signal, i);
				minimum = minimum + ( ( LinearValue - HistogramValue ) * ( LinearValue - HistogramValue ) );
			}

			for(int i=index_M; i < Histogram_size; i++)
			{
				auto HistogramValue = gsl_vector_int_get(Histogram->signal, i);
				minimum = minimum + (HistogramValue*HistogramValue);
			}

			if(minimum < globalminimum)
			{
				globalminimum=minimum;
				N=index_N;
				M=index_M;
			}
			minimum=0;
		}
	}

	gsl_vector_int_memcpy(this->histogram_y->signal,Histogram->signal);
	this->N=N;
	this->M=M;
	this->X=X;
	this->Y=Y;
	this->HRVTriangularIndex = RR_intervals_size/Y;
	this->TINN = M - N;

}

void GeometricAnalysis::MakePoincareAndSDParams()
{
	IntSignal poincare_x;
	IntSignal poincare_y;
	IntSignal diff;

	auto RR_intervals_size = RR_intervals->signal->size;

	double *data_SDNN = new double[RR_intervals_size];

	for(int i = 0; i < RR_intervals_size; i++) 
	{
		auto SignalValue = gsl_vector_get (RR_intervals->signal, i);
		data_SDNN[i] = SignalValue;

		if(i < RR_intervals_size-1) gsl_vector_int_set(poincare_x->signal, i, SignalValue);
		if(i > 0) gsl_vector_int_set(poincare_y->signal, i, SignalValue);
	}

	double sdnn = gsl_stats_sd(data_SDNN,1,RR_intervals_size);

	gsl_vector_int_memcpy(diff->signal, this->poincare_x->signal);
	gsl_vector_int_sub(diff->signal,poincare_y->signal);

	double *data_DIFF = new double[RR_intervals_size];

	for(int i = 0; i < RR_intervals_size; i++) 
	{
		auto SignalValue = gsl_vector_int_get (diff->signal, i);
		data_DIFF[i] = SignalValue;
	}

	double sdsd = gsl_stats_sd(data_DIFF,1,RR_intervals_size);

	double SD1 = sqrt(0.5*pow(sdsd,2));
	double SD2 = sqrt(2*pow(sdnn,2) - 0.5*pow(sdsd,2));
	this->SD1 = SD1;
	this->SD2 = SD2;
	gsl_vector_int_memcpy(this->poincare_x->signal,poincare_x->signal);
	gsl_vector_int_memcpy(this->poincare_y->signal,poincare_y->signal);
}

void GeometricAnalysis::SetHRV2Params()
{
	hrv2.SetSD1(this->SD1);
	hrv2.SetSD2(this->SD2);
	hrv2.SetTINN(this->TINN);
	hrv2.SetM(this->M);
	hrv2.SetN(this->N);
	hrv2.SetHRVTriangularIndex(this->HRVTriangularIndex);
	hrv2.SetY(this->Y);
	hrv2.SetX(this->X);
	hrv2.SetHistogramBinLength(this->HistogramBinLength);

}