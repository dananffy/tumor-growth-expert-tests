function result = T01_demo(mode)
% T01_DEMO  在 MATLAB 中运行 T01 最小实验、计算参考解并显示结果。
%   T01_demo             新计算一次并绘图（默认；可直接点击编辑器“运行”）。
%   T01_demo('prepare')  会前预编译；不运行数值实验。
%   T01_demo('replay')   重开最近一次结果，明确标注为已有结果回放。
%
% MATLAB: 控制、数据检查、独立 ode45 参考解、绘图和保存。
% 二维自由边界求解器: 原已核验 C++ 程序，经现有 WSL Ubuntu-24.04 调用。
% Windows 端无需 Python 科学计算包；不需要额外 MATLAB 工具箱。
% 固定输入: N=64, dt=0.04, T=0.2, 32 控制点；用于现场流程与径向参考演示。

if nargin == 0, mode = 'run'; end
mode = validatestring(mode, {'run','prepare','replay'});
demoRoot = fileparts(mfilename('fullpath'));
outputRoot = fullfile(demoRoot, 'runs');
if ~isfolder(outputRoot), mkdir(outputRoot); end
if strcmp(mode, 'replay')
    pointer = fullfile(outputRoot, 'latest.txt');
    assert(isfile(pointer), 'T01:NoResult', '尚无结果，请先运行 T01_demo。');
    outputDir = fullfile(outputRoot, strtrim(fileread(pointer)));
    loaded = load(fullfile(outputDir,'T01_result.mat'), 'data','summary');
    fig = draw_result(loaded.data, loaded.summary, true);
    result = struct('output_directory',outputDir,'summary',loaded.summary,'figure',fig);
    fprintf('已有结果回放（没有新计算）：\n%s\n', outputDir);
    return
end
assert(ispc, 'T01:Platform', '本入口为本机 Windows + WSL 环境准备。');
backend = fullfile(demoRoot, 'backend');
driver = fullfile(backend, 'frozen_cpp_runner.py');
driverWSL = to_wsl_path(driver);
command = sprintf('wsl.exe -d Ubuntu-24.04 -- python3 "%s" --case T01_radial_64', driverWSL);
if strcmp(mode,'prepare'), command = [command ' --build-only']; end
started = tic;
fprintf('\nT01 圆形肿瘤径向生长：%s\n', mode);
fprintf('固定配置：64×64 网格，32 控制点，dt=0.04，T=0.2\n');
fprintf('lambda=0.5，c_B=10，G0=0.2，R0=0.5\n');
fprintf('正在调用求解器；首次运行需编译，后续自动复用编译缓存。\n');
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
outputDir = fullfile(outputRoot, [stamp '_' mode]);
mkdir(outputDir);
[status, consoleText] = system(command, '-echo');
write_utf8(fullfile(outputDir, 'execution.log'), consoleText);
if status ~= 0
    error('T01:SolverFailed', '求解器退出码 %d；日志：%s', status, fullfile(outputDir,'execution.log'));
end
if strcmp(mode,'prepare')
    result = struct('mode','prepare','elapsed_seconds',toc(started),'output_directory',outputDir);
    write_utf8(fullfile(outputDir,'prepare.json'), jsonencode(result,PrettyPrint=true));
    fprintf('会前准备完成，用时 %.2f 秒。现场直接运行 T01_demo。\n',result.elapsed_seconds);
    return
end

runID = strtrim(fileread(fullfile(backend,'runs','T01','latest_run.txt')));
rawDir = fullfile(backend,'runs','T01',runID);
meta = jsondecode(fileread(fullfile(rawDir,'config_and_environment.json')));
assert(meta.exit_code == 0 && strcmp(meta.case_id,'T01_radial_64'), 'T01:Metadata', '求解记录不匹配。');
assert(strcmp(runID,meta.run_id), 'T01:RunID', '运行目录与运行记录不匹配。');
raw = readtable(fullfile(rawDir,'radius_diagnostics.csv'),VariableNamingRule='preserve');
assert(all(isfinite(raw{:,:}), 'all'), 'T01:Nonfinite', '原始半径记录存在非有限数。');
t = raw.time;
expected = (0:meta.config.dt:meta.config.final_time)';
assert(numel(t)==numel(expected) && max(abs(t-expected))<1e-12, 'T01:Times', '保存时刻不符合固定配置。');
params = meta.parameters;
rhs = @(~,r) params.G0*params.c_B .* besseli(1,sqrt(params.lambda)*r) ./ ...
    (sqrt(params.lambda)*besseli(0,sqrt(params.lambda)*r));
% 此处独立用 MATLAB 积分，不读入预先保存的参考曲线。
reference = ode45(rhs, [0 meta.config.final_time], params.R0, ...
    odeset(RelTol=1e-11,AbsTol=1e-13));
tight = ode45(rhs, [0 meta.config.final_time], params.R0, ...
    odeset(RelTol=1e-12,AbsTol=1e-14));
Rref = deval(reference,t)';
Rcheck = deval(tight,t)';
assert(all(isfinite(Rref)) && all(Rref>0), 'T01:Reference', '参考积分失败。');
points = cell(numel(t),1);
for k=1:numel(t)
    points{k} = readmatrix(fullfile(rawDir,'result','points',sprintf('ctrl_points_%d.txt',k-1)));
    assert(isequal(size(points{k}),[params.n_ctrl,2]) && all(isfinite(points{k}),'all'), ...
        'T01:Points','边界控制点文件不完整。');
    for field = {'c','p'}
        values = readmatrix(fullfile(rawDir,'result','functions',sprintf('%s_%d.txt',field{1},k-1)));
        assert(isequal(size(values),[meta.config.grid+1,meta.config.grid+1]) && all(isfinite(values),'all'), ...
            'T01:Field','营养或压力数据不完整。');
    end
end
logText = fileread(fullfile(rawDir,'run.log'));
assert(~any(contains(logText,{'failed to converge','Invalid velocity','Could not open file'})), ...
    'T01:SolverLog','求解日志存在失败信息。');
absoluteError = abs(raw.radius_spline-Rref);
relativeError = absoluteError./Rref;
comparison = table(t,raw.radius_spline,Rref,absoluteError,relativeError, ...
    VariableNames={'time','radius_numerical','radius_reference_matlab','absolute_error','relative_error'});
writetable(comparison,fullfile(outputDir,'radius_comparison_matlab.csv'));
data = struct('t',t,'radius',raw.radius_spline,'reference',Rref, ...
    'absolute_error',absoluteError,'relative_error',relativeError,'points',{points}, ...
    'dense_t',linspace(0,t(end),300),'params',params,'grid',meta.config.grid,'dt',meta.config.dt);
data.dense_reference = deval(reference,data.dense_t);
data.pressure = readmatrix(fullfile(rawDir,'result','functions',sprintf('p_%d.txt',numel(t)-1)));
data.nutrient = readmatrix(fullfile(rawDir,'result','functions',sprintf('c_%d.txt',numel(t)-1)));
data.coords = linspace(params.box(1),params.box(2),meta.config.grid+1);
summary = struct('run_id',runID,'matlab_release',version('-release'), ...
    'recorded_at_local',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'mode','new_solver_run_and_matlab_reference','grid',meta.config.grid,'dt',meta.config.dt, ...
    'final_time',t(end),'saved_states',numel(t),'radius_numerical_final',raw.radius_spline(end), ...
    'radius_reference_final',Rref(end),'absolute_radius_error_final',absoluteError(end), ...
    'relative_radius_error_final',relativeError(end), ...
    'reference_tolerance_crosscheck',max(abs(Rref-Rcheck)), ...
    'all_checked_numeric_arrays_finite',true,'solver_seconds',meta.run_seconds, ...
    'build_seconds',meta.build_seconds,'compute_and_reference_seconds',toc(started), ...
    'source_metadata',meta,'scope','T01 short radial demonstration; not a full refinement test or expert acceptance conclusion.');
% 独立保留这一轮原始数据，不依赖 latest 指针查验历史结果。
copyfile(rawDir,fullfile(outputDir,'raw'));
write_utf8(fullfile(outputDir,'summary.json'),jsonencode(summary,PrettyPrint=true));
save(fullfile(outputDir,'T01_result.mat'),'data','summary');
write_utf8(fullfile(outputRoot,'latest.txt'),[stamp '_' mode]);
fig = draw_result(data,summary,false);
exportgraphics(fig,fullfile(outputDir,'T01_matlab_demo.png'),Resolution=160);
exportgraphics(fig,fullfile(outputDir,'T01_matlab_demo.pdf'),ContentType='vector');
savefig(fig,fullfile(outputDir,'T01_matlab_demo.fig'));
summary.total_seconds = toc(started);
write_utf8(fullfile(outputDir,'summary.json'),jsonencode(summary,PrettyPrint=true));
save(fullfile(outputDir,'T01_result.mat'),'data','summary');
fprintf('\n新计算完成：%s\n',runID);
fprintf('终点数值半径 %.10f，MATLAB参考半径 %.10f\n',raw.radius_spline(end),Rref(end));
fprintf('终点相对误差 %.4f%%；参考积分容差复核差 %.3g\n',100*relativeError(end),summary.reference_tolerance_crosscheck);
fprintf('结果已保存：\n%s\n',outputDir);
result = struct('output_directory',outputDir,'summary',summary,'figure',fig);
end

function fig = draw_result(d,s,replay)
if replay, label='已有结果回放（本次未计算）'; else, label='本次新计算'; end
fig = figure(Name=['T01 MATLAB 演示 — ' label],Color='w',Position=[80 80 1280 820]);
layout = tiledlayout(fig,2,2,TileSpacing='loose',Padding='loose');
theta = linspace(0,2*pi,300);
ax1 = nexttile(layout);
hold(ax1,'on');
plot(ax1,d.params.R0*cos(theta),d.params.R0*sin(theta),':',Color=[.55 .55 .55],DisplayName='初始边界');
plot(ax1,d.reference(end)*cos(theta),d.reference(end)*sin(theta),'k-',LineWidth=1.7,DisplayName='终点参考圆');
q = d.points{end}; q=[q;q(1,:)];
plot(ax1,q(:,1),q(:,2),'o-',Color=[.1 .4 .75],MarkerSize=4,DisplayName='数值边界控制点');
axis(ax1,'equal'); grid(ax1,'on'); xlabel(ax1,'x'); ylabel(ax1,'y');
title(ax1,sprintf('边界推进：t=0 → %.2f',d.t(end))); legend(ax1,Location='best');
ax2 = nexttile(layout);
plot(ax2,d.dense_t,d.dense_reference,'k-',LineWidth=1.6,DisplayName='MATLAB ode45 参考'); hold(ax2,'on');
plot(ax2,d.t,d.radius,'o--',Color=[.1 .4 .75],LineWidth=1.3,DisplayName='二维求解器半径');
grid(ax2,'on'); xlabel(ax2,'时间 t'); ylabel(ax2,'半径 R(t)'); title(ax2,'半径随时间的变化'); legend(ax2,Location='best');
ax3 = nexttile(layout);
plot(ax3,d.t,100*d.relative_error,'o-',Color=[.8 .25 .12],LineWidth=1.5);
grid(ax3,'on'); xlabel(ax3,'时间 t'); ylabel(ax3,'相对误差（%）');
title(ax3,sprintf('终点半径相对误差：%.4f%%',100*s.relative_radius_error_final));
ax4 = nexttile(layout);
[xx,yy]=ndgrid(d.coords,d.coords);
pressure=d.pressure;
radialMin=min(vecnorm(d.points{end},2,2));
pressure(hypot(xx,yy)>radialMin)=NaN;
h=imagesc(ax4,d.coords,d.coords,pressure'); set(h,AlphaData=isfinite(pressure'));
set(ax4,YDir='normal'); axis(ax4,'image'); hold(ax4,'on');
plot(ax4,q(:,1),q(:,2),'k-',LineWidth=1);
colorbar(ax4); xlim(ax4,[-.75 .75]); ylim(ax4,[-.75 .75]); xlabel(ax4,'x'); ylabel(ax4,'y'); title(ax4,'终点压力场 p(x,y)');
title(layout,sprintf('%s｜T01 径向生长｜N=%d，dt=%.2f，T=%.2f',label,d.grid,d.dt,d.t(end)),FontWeight='bold');
subtitle(layout,'MATLAB 控制、独立参考与绘图；C++/WSL 二维求解。');
set(findall(fig,'Type','axes'),FontSize=12,FontName='Microsoft YaHei');
set(findall(fig,'Type','text'),FontName='Microsoft YaHei');
drawnow;
end

function p = to_wsl_path(p)
p = strrep(p,'\','/');
assert(numel(p)>2 && p(2)==':','T01:Path','演示目录应位于本地盘符路径。');
p = ['/mnt/' lower(p(1)) p(3:end)];
end

function write_utf8(file,text)
fid=fopen(file,'w','n','UTF-8');
assert(fid>=0,'T01:Write','无法写入文件：%s',file);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',text);
end
