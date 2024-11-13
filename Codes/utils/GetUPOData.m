function [xdat,PeriodData] = GetUPOData(folder, nUPO,Case,t,dt,ode_options,LORENZ)
x0 = zeros(3,nUPO);
k=0;
PeriodData = zeros(nUPO,2);
if(strcmp(Case,'Case1'))
    for counter=1:nUPO/2
        k=k+1;
        str = sprintf('orbit%d.dat',counter);
        filename = append(folder,str);
        orbit=load(filename);
        x0 = orbit(:,1);
        [~,xdat(:,:,k)] = ode45(@(t,x)LORENZ(t,x),t,x0,ode_options);
        str = sprintf('T%d.dat',counter);
        filename = append(folder,str);
        PeriodData(k,1) =  load(filename);
        PeriodData(k,2) = ceil(PeriodData(k,1)/dt);
        k=k+1;
        x0(1) = x0(1)*-1;
        x0(2) = x0(2)*-1;
        [~,xdat(:,:,k)] = ode45(@(t,x)LORENZ(t,x),t,x0,ode_options);
        PeriodData(k,1) = load(filename);
        PeriodData(k,2) = ceil(PeriodData(k,1)/dt);
        disp(counter);
    end
elseif(strcmp(Case,'Case2')||strcmp(Case,'Case3'))
    for counter=1:nUPO
        k=k+1;
        str = sprintf('orbit%d.dat',counter);
        filename = append(folder,str);
        orbit=load(filename);
        x0 = orbit(:,1);
        [~,xdat(:,:,k)] = ode45(@(t,x)LORENZ(t,x),t,x0,ode_options);
        str = sprintf('T%d.dat',counter);
        filename = append(folder,str);
        PeriodData(k,1) =  load(filename);
        PeriodData(k,2) = ceil(PeriodData(k,1)/dt);
        disp(counter);
    end

end
