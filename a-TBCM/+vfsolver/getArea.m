function [Amin,Pmin] = getArea(x,model,Lg,x_du0,x_dl0)
  switch upper(model)
    case {'BCM','BCM+'}     
      Amin = 2*Lg*max(0,min(x(1),x(2)));
      Pmin = 2*(Lg*double(max(0,min(x(1),x(2)))>0) + max(0,min(x(1),x(2))));
    case 'TBCM'
      % Original code (time consuming)
      %   Border = max(0,min(x(1)+linspace(0,x_du0,1e5),x(2)+linspace(0,x_dl0,1e5)));
      %   Amin = 2*Lg*sum(Border)/1e5;
      %   Pmin = 2*sum(sqrt(diff(Border).^2+(Lg/1e5)^2).*double(Border(2:end)>0))+Border(1)+Border(end);
      % New Code
      Border(1) = min(1,max(0,-x(1)/x_du0)); % Contact of the upper mass (normalized by Lg)
      Border(2) = min(1,max(0,-x(2)/x_dl0)); % Contact of the lower mass
      Border(3) = min(1,max(0,(x(1)-x(2))/(x_dl0-x_du0))); % Contact between masses
      Border(4) = 0;
      Border(5) = 1;

      Border = sort(Border,2,'ascend');
      % Displacements in the critical points
%       Border(2,:) = [...
%         max(0,min(x(1)+x_du0*Border(1),x(2)+x_dl0*Border(1)));...
%         max(0,min(x(1)+x_du0*Border(2),x(2)+x_dl0*Border(2)));...
%         max(0,min(x(1)+x_du0*Border(3),x(2)+x_dl0*Border(3)));...
%         max(0,min(x(1)+x_du0*Border(4),x(2)+x_dl0*Border(4)));...
%         max(0,min(x(1)+x_du0*Border(5),x(2)+x_dl0*Border(5)))];
      
      Border(2,:) = max(0,min(x(1) + x_du0*Border,x(2) + x_dl0*Border));

      % Areas of the different sections
%       Areas = [...
%         (Border(1,2) - Border(1,1))*(Border(2,1) + Border(2,2))/2;...
%         (Border(1,3) - Border(1,2))*(Border(2,2) + Border(2,3))/2;...
%         (Border(1,4) - Border(1,3))*(Border(2,3) + Border(2,4))/2;...
%         (Border(1,5) - Border(1,4))*(Border(2,4) + Border(2,5))/2];
      Area = (Border(1,2:5) - Border(1,1:4)).*(Border(2,2:5) + Border(2,1:4))/2;

      % Perimeters of the differet sections
%        Perimeter = [...
%         Border(2,1);...
%         sqrt(((Border(1,2) - Border(1,1))*Lg)^2 + (Border(2,2)-Border(2,1))^2)*(Areas(1)>0);...
%         sqrt(((Border(1,3) - Border(1,2))*Lg)^2 + (Border(2,3)-Border(2,2))^2)*(Areas(2)>0);...
%         sqrt(((Border(1,4) - Border(1,3))*Lg)^2 + (Border(2,4)-Border(2,3))^2)*(Areas(3)>0);...
%         sqrt(((Border(1,5) - Border(1,4))*Lg)^2 + (Border(2,5)-Border(2,4))^2)*(Areas(4)>0);...
%         Border(2,5)];
      Perimeter = [Border(2,1) sqrt((Border(1,2:5) - Border(1,1:4)).^2 + (Border(2,2:5) - Border(2,1:4)).^2).*(Area>0) Border(2,5)];

      Amin = 2*sum(Area)*Lg;
      Pmin = 2*sum(Perimeter);
  end
end