function [fn,psd,z] = freevib(x,varargin)
%{
DESCRIPTION
This function estimates the vibration properties of tree parts using
measurements of tree movement undergoing free vibration.

NOTES
The function accepts tree displacement and acceleration time histories 
recorded in SI units. A multiple or submultiple of the SI unit can be 
formed by specifying the appropriate SI unit prefix. Free vibration 
time histories must be trimmed to start at a local minimum or maximum.
Damping ratio can be estimated using the log-decrement method ('decrement')
or by fitting the general solution to the equation of motion for a damped
harmonic oscillator ('curve'). 

[fn,psd,z] = freevib(x);
[fn,psd,z] = freevib(x,'parameter1',value1,...);

REQUIRED INPUTS
x: float - m x n matrix with tree response measurements

OPTIONAL INPUTS
fsample: float - sampling frequency (Hz)
t: float - elapsed time (sec) (for uneven sample intervals)
quantity: string - measurement quantity, i.e., 'displacement' or
'acceleration'
siprefix: string - SI unit prefix, e.g., 'deci', 'milli', 'micro'
damping: string - 'decrement' or 'curve'

PARAMETER                     CLASS       DEFAULT VALUE
---------------------------------------------------------------------------
fsample                       float       1
t                             float       (0:length(x)-1)'
quantity                      string      'displacement'
siprefix                      string      ''
damping                       string      'decrement'

EXAMPLE
t=(0:0.05:60)';
y=10.*exp(-0.05*(0.2*2*pi/(sqrt(1-0.05.^2))*t)).*sin(0.2*2*pi.*t+(pi/2));

OUTPUTS
fn: Natural frequency
psd: Power spectral density
z: Damping ratio
%}
format long

%SI unit prefixes
A={'giga', 'G'; 'mega', 'M'; 'kilo', 'k'; 'hecto', 'h'; 'deka', 'da';...
    'deci', 'd'; 'centi', 'c'; 'milli', 'm'; 'micro', '\mu'; 'nano', 'n'};

% Check function call
if nargin < 1
    error('treevib requires at least one input');
elseif ~isempty(varargin)
    if mod(length(varargin),2)~=0
        if (length(varargin)==1 && isstruct(varargin{1}))
            varargin = reshape([fieldnames(varargin{1})...
                struct2cell(varargin{1})]',1,[]);
        else
            error(strcat('Inputs must be paired: treevib(x,',...
                '''PropertyName'',PropertyValue,...)'));
        end
    elseif any(ismember(varargin(1:2:size(varargin,2)),'damping'))
        if ~any(strcmp(varargin{find(strcmp(varargin,'damping'))+1},...
            {'decrement','curve'}))
                error('damping must be set equal to ''decrement'' or ''curve''');
        end
    end
end

% Default parameters
fs = 1;
s = (0:length(x)-1)';
response = 'displacement';
prefix = '';
zeta = 'decrement';

% User-specified parameters
if ~isempty(varargin)
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'fsample'
                fs = varargin{i+1};
            case 't'
                s = varargin{i+1};
            case 'quantity'
                response = varargin{i+1};
            case 'siprefix'
                prefix = varargin{i+1};
            case 'damping'
                zeta = varargin{i+1};
        otherwise
            error(strcat([varargin{i} ' is not a valid property for',...
                'the treevib function.']));
        end
    end
end

%If x is a matrix, make a projection onto the resultant vector
if size(x,2)>1
    x=dot(x,repmat(mean(x),length(x),1),2)./norm(mean(x));
end

%Frequency estimate
%Detrend time history and remove mean
x=x-mean(x);
x=detrend(x);

%Compute spectral estimate
if any(ismember(varargin(1:2:size(varargin,2)),'t'))
    if abs(mean(diff(diff(s)))) < 1e-10 %Welch's PSD Estimate
        dt=mean(diff(s));
        fs=1/dt;
        [pxx,f]=pwelch(x,length(x),0,length(x),fs);
    elseif abs(mean(diff(diff(s)))) >= 1e-10 %Lomb-Scargle Periodogram
        [pxx,f]=plomb(x,s);
    end
elseif ~any(ismember(varargin(1:2:size(varargin,2)),'t'))
    s=(0:1/fs:(length(x)-1)*1/fs)';
    [pxx,f]=pwelch(x,length(x),0,length(x),fs);
end

%Find peak frequency(ies)
[G,H,~,I]=findpeaks(pxx);
J=[G,H,I];
[~,K]=sort(J(:,3),'descend');
J=J(K,:);
L=J(:,3)>=(0.3*J(1,3));
J=J(L,:);
fd=f(J(:,2));
psd=J(:,1);

%Damping ratio estimate
z=zeros(size(fd));
yr=zeros(length(x),size(fd,1));

for i=1:size(fd,1)
    if strcmp(zeta,'decrement')
        %Logarithmic decrement
        [pks,~]=findpeaks(x,s,'MinPeakDistance',0.85*1/fd(i));
        delta=polyfit((1:numel(pks))',log(pks),1);
        z(i)=abs(delta(1)/sqrt(4*pi.^2+delta(1).^2));
    elseif strcmp(zeta,'curve')
        %Damped sine curve fit
        [~,~,~,z(i),~,yr(:,i)]=sf_engine(x,s,length(x),fd(i)-0.01,fd(i)+0.01,...
            1e6,s(end)-s(1),0,0,0,0,0);
    end
end

%Determine natural frequency
fn=fd./sqrt(1-z.^2);

%Adjust plot y-axis
M=abs(max(x));
ymax=(roundn(2*M,-1)*0.5)+0.05;
ymin=-1*((roundn(2*M,-1)*0.5)+0.05);
N=roundn(ymax/2,-1);
O=N*[-1 0 1];
P=1.1*max(pxx);

% Create figure
figure1 = figure('Color',[1 1 1],'Units','centimeters',...
    'Position',[1 1 27 14],'PaperPositionMode','Auto');

% Create subplot
subplot1 = subplot(2,1,1,'Parent',figure1,'XGrid','on','XScale','log');
xlim(subplot1,[0.1 10]);
ylim(subplot1,[0 P]);
box(subplot1,'on');
hold(subplot1,'on');

% Create semilogx
semilogx(f,pxx,'Parent',subplot1,'LineWidth',1,'Color',[0 0 0]);

% Create text
for i=1:numel(fd)
    text('Parent',subplot1,'FontWeight','bold','FontSize',13,'FontName',...
        'Helvetica','String',[' \leftarrow ',num2str(roundn(fd(i),-2)),' '],...
        'Position',[fd(i)+0.005 psd(i) 0]);
end

% Create xlabel
xlabel('Frequency (Hz)','FontSize',13);

% Create ylabel
if strcmp(response,'displacement')
    ylabel(strcat('Power/Frequency (',A{ismember(A(:,1),prefix,'legacy'),2},...
        'm^{2}\cdotHz^{-1})'),'FontSize',13);
elseif strcmp(response,'acceleration')
    ylabel(strcat('Power/Frequency (',A{ismember(A(:,1),prefix,'legacy'),2},...
        'm^{2}\cdots^{-4}\cdotHz^{-1})'),'FontSize',13);
end

% Create title
title('Power Spectral Density','Units','normalized',...
    'HorizontalAlignment','right',...
    'FontWeight','bold',...
    'FontSize',13,...
    'VerticalAlignment','middle',...
    'Position',[1 1.06 0]);

% Create subplot
subplot2 = subplot(2,1,2,'Parent',figure1,'YGrid','on','XGrid','on',...
    'YLim',[ymin ymax],'YTickLabel',O,'YTick',O);
box(subplot2,'on');
hold(subplot2,'on');

% Create plot
plot(s,x,'Parent',subplot2,'LineWidth',0.5,'Color',[0 0 0]);
hold on;
if strcmp(zeta,'decrement')
    for i=1:numel(fd)
        plot(s,abs(x(1)).*exp(-z(i).*(fn(i)*2*pi).*s),'Parent',subplot2,...
            'Linewidth',2);
    end
elseif strcmp(zeta,'curve')
    for i=1:numel(fd)
        plot(s,yr(:,i),'Parent',subplot2,'Linewidth',2);
    end
end
for i=1:numel(fd)
    text('Parent',subplot2,'FontWeight','bold','FontSize',13,'FontName',...
        'Helvetica','String',['\zeta = ',num2str(roundn(z(i),-2)),' '],...
        'Position',[0.9*s(end) 0.85*ymax 0]);
end
hold off;

% Create xlabel
xlabel('Time (s)','FontSize',13);

% Create ylabel
if strcmp(response,'displacement')
    ylabel(strcat('Amplitude (',A{ismember(A(:,1),prefix,'legacy'),2}...
        ,'m)'),'FontSize',13);
elseif strcmp(response,'acceleration')
    ylabel(strcat('Amplitude (',A{ismember(A(:,1),prefix,'legacy'),2}...
        ,'m\cdots^{-2})'),'FontSize',13);
end

% Create title
title(strcat(upper(response(1)),response(2:end),' Time History'),...
    'Units','normalized','HorizontalAlignment','right',...
    'FontWeight','bold','FontSize',13,'VerticalAlignment','middle',...
    'Position',[1 1.06 0]);
end

function [x1r,x2r,x3r,x4r,x5r,yr] = sf_engine(a,t,num2,flow,fup,nt,dur,...
    x1r,x2r,x3r,x4r,x5r)

rng(10,'v4');
tp=2*pi;

ave=mean(a);
sd=std(a);

am=2.*sd;
n=num2;
pu=1.;
pl=0.;

errormax=1.0e+53;

jk=0;

delta=0.001;
nnn=fix(0.1*nt);

for j=1:nt

    if(j>nnn)
      fr=(x2r/tp);
      if( (fup-fr)/fup < delta)
          fup=fup*(1+delta);
      end
      if(flow>1.0e-06)
        if( (fr-flow)/flow < delta)
            flow=flow*(1-delta);
        end  
      end
    end

	jk=jk+1;
      
    if(jk==10000)
		jk=0;
    end

	x1=rand;
	x2=rand;
	x3=rand;
	x4=rand;
	x5=rand;

    n1=fix(2.*nt/3.);
    n2=fix(4.*nt/5.);
    if( j < n1 )
		x1=2.*am*x1;
		x2=((fup-flow)*x2+flow)*tp;
		x3=((pu-pl)*x3+pl)*tp;
		x4=(x4^2.5);
		x5=dur*(x5^2.);
    end
    if( j >= n1 && j < n2)
		x1=x1r*(0.98+0.04*x1);
		x2=x2r*(0.98+0.04*x2);
		x3=x3r*(0.98+0.04*x3);
		x5=x5r*(0.98+0.04*x5);
    end
    if( j >= n2 )
		x1=x1r*(0.99+0.02*x1);
		x2=x2r*(0.99+0.02*x2);
		x3=x3r*(0.99+0.02*x3);
		x4=x4r*(0.90+0.20*x4);
		x5=x5r*(0.99+0.02*x5);
    end

    if x2 > fup*tp
        x2 = fup*tp;
    end

    if x2 < flow*tp
        x2 = flow*tp;
    end    

    if x4 >= 0.85
        x4=0.05;
    end    

    if j == 1
        x1=0.;
    elseif x1 < 1.0e-12
        x1=am*rand;
    end

	error=0.;
       
    domegan=x4*x2;
    omegad=x2*sqrt(1.-x4^2);
	
    yr=zeros(n,1);
    for i=1:n
		tt=t(i)-t(1);
		y=0.;
        if ( tt > x5 )
            ta=tt-x5;
            y=x1*exp(-domegan*ta)*sin(omegad*ta-x3);
            yr(i)=y;
        end	
		error=error+((a(i)-y)^2.);
    end
	error=sqrt(error);
    if ( error < errormax )
		x1r=x1;
		x2r=x2;
		x3r=x3;
		x4r=x4;
		x5r=x5;
	    errormax=error;
    end				
end

ave=0.;
sq=0.;

for i=1:n
    tt=t(i)-t(1);
    domegan=x4r*x2r;
    omegad=x2r*sqrt(1.-x4r^2);
    if(tt>x5r)
        ttt=tt-x5r;
        a(i)=a(i)-x1r*exp( -domegan*ttt )*sin( omegad*ttt - x3r );
        ave=ave+a(i);
	    sq=sq+(a(i)^2.);
    end
end

end