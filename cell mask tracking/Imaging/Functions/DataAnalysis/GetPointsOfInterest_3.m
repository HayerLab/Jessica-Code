function [POI,badtraces]=GetPointsOfInterest_3(signal,traces,stats,minlength)
switch signal
    case 'Cdk2'
        %[POI,badtraces]=getCdk2features_1(traces,stats,minlength);
        %[POI,badtraces]=getCdk2features_steve(traces,stats,minlength);
        [POI,badtraces]=getCdk2features_2(traces,stats,minlength);
    case 'Cdk2_SR'
        [POI,badtraces]=getCdk2features_SR(traces,stats,minlength);
    case 'Geminin'
        [POI,badtraces]=getGemininfeatures_1(traces,stats);
    case 'Geminin_SR'
        [POI,badtraces]=getGemininfeatures_serumrelease(traces,stats);
    case 'PipdegFall'
        %[POI,badtraces]=getPipdegFall_NormByFirst(traces,stats,minlength);
        %[POI,badtraces]=getPipdegFall_NormByMother(traces,stats,minlength);
        [POI,badtraces]=getPipdegFall_Dropstart(traces,stats,minlength);
    case 'PipdegRise'
        [POI,badtraces]=getPipdegRise(traces,stats,minlength);
end