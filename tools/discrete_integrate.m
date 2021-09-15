function integral = discrete_integrate(xData, yData)

%     function integral = discrete_integrate(xData, yData)
%
%           dXData = diff(xData);
%           dIntegral = dXData.*(yData(1:end-1)+yData(2:end))/2;
%           integral = sum(dIntegral);

    
    dXData = diff(xData);
    dIntegral = dXData.*(yData(1:end-1)+yData(2:end))/2;
    integral = sum(dIntegral);