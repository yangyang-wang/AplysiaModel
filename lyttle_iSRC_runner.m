% Generate ISRC with global uniform or piecewise uniform rescaling; start from close
pt=true;
ICisClose=true;
isUniform=true;

% iSRC with uniform rescaling, starting from close
[t_isrc_uni_close,x_isrc_uni_close]=lyttle_iSRC(pt,ICisClose,isUniform);

% iSRC with uniform rescaling, starting from open
[t_isrc_uni_open,x_isrc_uni_open]=lyttle_iSRC(pt,~ICisClose,isUniform);

% iSRC with piecewise rescaling, starting from close
[t_isrc_pw_close,x_isrc_pw_close]=lyttle_iSRC(pt,ICisClose,~isUniform);

% iSRC with piecewise rescaling, starting from open
[t_isrc_pw_open,x_isrc_pw_open]=lyttle_iSRC(pt,~ICisClose,~isUniform);

% text(-1,1.05,'$\rm (A)$','Interpreter','latex','FontSize',20,'FontWeight','bold')
