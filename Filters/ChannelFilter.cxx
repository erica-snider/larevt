///////////////////////////////////////////////////////
//
// ChannelFilter Class
// 
// This class has been obsoleted and is now a deprecated interface for
// IChannelFilterService.
// 
// Please update your code to use the service directly.
// 
// 
// Original class: pagebri3@msu.edu
//
///////////////////////////////////////////////////////


// Our header
#include "Filters/ChannelFilter.h"


// Framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Utilities/Exception.h"

// LArSoft libraries
#include "CalibrationDBI/Interface/IChannelFilterService.h"
#include "CalibrationDBI/Interface/IChannelFilterProvider.h"


filter::ChannelFilter::ChannelFilter() {
  
  if (!&*(art::ServiceHandle<lariov::IChannelFilterService>())) {
    throw art::Exception(art::errors::Configuration)
      << "Failed to obtain an instance of IChannelFilterService service"
      ;
  }
  LOG_ERROR("ChannelFilter") << "ChannelFilter is now deprecated."
    " Replace it with IChannelFilterService";
  
} // filter::ChannelFilter::ChannelFilter()


///////////////////////////////////////////////////////
bool filter::ChannelFilter::BadChannel(uint32_t channel) const {
  return art::ServiceHandle<lariov::IChannelFilterService>()
    ->GetFilter().IsBad(channel);
}

///////////////////////////////////////////////////////
bool filter::ChannelFilter::NoisyChannel(uint32_t channel) const{
  return art::ServiceHandle<lariov::IChannelFilterService>()
    ->GetFilter().IsNoisy(channel);
}

///////////////////////////////////////////////////////
std::set<uint32_t> filter::ChannelFilter::SetOfBadChannels() const {
  return art::ServiceHandle<lariov::IChannelFilterService>()
    ->GetFilter().BadChannels();
}

///////////////////////////////////////////////////////
std::set<uint32_t> filter::ChannelFilter::SetOfNoisyChannels() const {
  return art::ServiceHandle<lariov::IChannelFilterService>()
    ->GetFilter().NoisyChannels();
}

///////////////////////////////////////////////////////
filter::ChannelFilter::ChannelStatus filter::ChannelFilter::GetChannelStatus(uint32_t channel) const
{
  
  lariov::IChannelFilterProvider const& filter
    = art::ServiceHandle<lariov::IChannelFilterService>()->GetFilter();
  
  if (!filter.IsPresent(channel)) return NOTPHYSICAL;
  if (filter.IsBad(channel)) return DEAD;
  if (filter.IsNoisy(channel)) return NOISY;
  return GOOD;
}
