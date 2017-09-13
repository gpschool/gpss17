%
% GPSS'17 appedix -
% Example of Kalman/RTS inference on GP regression problem
% with the Ornstein-Uhlenbeck covariance function. By SS'17.
%

    %%
    % Data
    %
    rng(20);
    
    X = 0:0.01:6;
    F = sin(X);
    
    s2 = 0.1^2;
    ind = sort(1+round(length(X) * rand(1,10)));

    x = X(ind)';
    y = F(ind)';
    y = y + sqrt(s2) * randn(size(y));
    
    clf;
    h = plot(X,F,x,y,'o');
    set(h,'LineWidth',2);
    
    %%
    % GP regression solution
    %
    q = 1;
    lam = 1;
    
    K_ou = @(x,xp) q/(2*lam) * exp(-lam * abs(x - xp));
    
    Kyy = zeros(length(x));
    Kfy = zeros(length(X),length(x));
    Kff = zeros(length(X),length(X));
    Sigma = s2 * eye(length(x));
    
    for i=1:size(Kyy,1)
        for j=1:size(Kyy,2)
            Kyy(i,j) = K_ou(x(i),x(j));
        end
    end
    Kyy = Kyy + Sigma;
    
    for i=1:size(Kfy,1)
        for j=1:size(Kfy,2)
            Kfy(i,j) = K_ou(X(i),x(j));
        end
    end
    
    for i=1:size(Kff,1)
        for j=1:size(Kff,2)
            Kff(i,j) = K_ou(X(i),X(j));
        end
    end
    
    mu = (Kfy / Kyy) * y;
    V  = Kff - (Kfy / Kyy) * Kfy';
    Vd = diag(V);
    
    QQ1 = mu-1.96*sqrt(Vd);
    QQ2 = mu+1.96*sqrt(Vd);

    clf;
    h2 = fill([X'; flipud(X')], [QQ1; flipud(QQ2)],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    hold on;
    
    h = plot(X,F,'k:',X,mu,x,y,'o');
    set(h,'LineWidth',2);
    
    %%
    % Kalman filter
    %
    MM = zeros(size(X));
    PP = zeros(size(X));

    m = 0;
    P = q / (2 * lam);
    
    Y = NaN * ones(size(X));
    Y(ind) = y;
    
    for k=1:length(X)
        if k > 1
            dt = X(k) - X(k-1);
            A = exp(-dt*lam);
            Q = q/(2*lam)*(1 - exp(-2*dt*lam));
            
            m = A * m;
            P = A^2 * P + Q;
        end
        
        if ~isnan(Y(k))
            K = P / (P + s2);
            m = m + K * (Y(k) - m);
            P = P - K^2 * (P + s2);

            QQ1 = MM-1.96*sqrt(PP);
            QQ2 = MM+1.96*sqrt(PP);

            clf;
            h2 = fill([X(1:k)'; flipud(X(1:k)')], [QQ1(1:k)'; flipud(QQ2(1:k)')],1);
            set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
            hold on;
            
            h = plot(X,F,'k:',X(1:k),MM(1:k),x,y,'o');
            set(h,'LineWidth',2);
            pause;
        end
        MM(k) = m;
        PP(k) = P;
    end
    
    QQ1 = MM-1.96*sqrt(PP);
    QQ2 = MM+1.96*sqrt(PP);

    clf;
    h2 = fill([X'; flipud(X')], [QQ1'; flipud(QQ2')],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])

    hold on;

    h = plot(X,F,'k:',X,MM,x,y,'o');
    set(h,'LineWidth',2);
    
    %%
    % RTS Smoother
    %
    MMS = MM;
    PPS = PP;
    
    ms = MM(end);
    Ps = PP(end);
    
    for k=length(MM)-1:-1:1
        dt = X(k+1) - X(k);
        A = exp(-dt*lam);
        Q = q/(2*lam)*(1 - exp(-2*dt*lam));

        m  = MM(k);
        P  = PP(k);
        mp = A * m;
        Pp = A^2 * P + Q;
        
        G = P * A' / Pp;
        
        ms = m + G * (ms - mp);
        Ps = P + G^2 * (Ps - Pp);
        
        MMS(k) = ms;
        PPS(k) = Ps;

        if ~isnan(Y(k))
            QQ1 = MMS-1.96*sqrt(PPS);
            QQ2 = MMS+1.96*sqrt(PPS);

            clf;
            h2 = fill([X'; flipud(X')], [QQ1'; flipud(QQ2')],1);
            set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])

            hold on;

            h = plot(X,F,'k:',X,MMS,x,y,'o');
            set(h,'LineWidth',2);
            
            pause;
        end
    end
    
    QQ1 = MMS-1.96*sqrt(PPS);
    QQ2 = MMS+1.96*sqrt(PPS);

    clf;
    h2 = fill([X'; flipud(X')], [QQ1'; flipud(QQ2')],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])

    hold on;

    h = plot(X,F,'k:',X,MMS,x,y,'o');
    set(h,'LineWidth',2);
            
    %%
    % Compare 1
    %
    clf;
    subplot(1,2,1);

    QQ1 = mu-1.96*sqrt(Vd);
    QQ2 = mu+1.96*sqrt(Vd);

    hold off;
    h2 = fill([X'; flipud(X')], [QQ1; flipud(QQ2)],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    hold on;
    
    h = plot(X,F,'k:',X,mu,x,y,'o');
    set(h,'LineWidth',2);
   
    
    subplot(1,2,2);
    hold off;

    QQ1 = MMS-1.96*sqrt(PPS);
    QQ2 = MMS+1.96*sqrt(PPS);

    h2 = fill([X'; flipud(X')], [QQ1'; flipud(QQ2')],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])

    hold on;

    h = plot(X,F,'k:',X,MMS,x,y,'o');
    set(h,'LineWidth',2);
    
    %%
    % Compare 2
    %
    clf;
    subplot(1,2,1);
    h = plot(X,mu,X,MMS,'--');
    set(h,'LineWidth',2);
    
    subplot(1,2,2);
    h = plot(X,Vd,X,PPS,'--');
    set(h,'LineWidth',2);

    max(abs(mu-MMS'))
    max(abs(Vd-PPS'))

    