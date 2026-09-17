from RayTracing.MacroConfig import RUNNING, ERROR, INIT, WARNING


def printMessage(case_name: str, status: str, expr: str):
    msg = '[' + case_name + ':' + status + ']  ' + expr
    print(msg)

    if status == ERROR:
        raise RuntimeError(expr)
    # elif status == WARNING:
    #     raise RuntimeWarning(expr)

def printLOGO():
    print('')
    print('')
    print('  ---------------------------------------------------------')
    print('        %%%%%%%%%  %%%%      %%%%%%%%%%%%%%     %%%%%%%%%%%')
    print('       %%     %%    %%%     %%%   %%     %%    %%   %%   %%')
    print('      %%      %%   %%%%   %% %    %%     %%    %%   %   %% ')
    print('      %%           % %%  %% %%   %%      %%   %%   %%   %% ')
    print('       %%%%       %% %% %%  %    %    %%%         %%       ')
    print('           %%%   %%  %%%%  %%   %%%%%%%           %%       ')
    print('            %%   %%  %%%  %%    %    %%          %%        ')
    print('   %%      %%   %%        %%   %%     %%         %%        ')
    print('  %%%     %%    %        %%   %%      %%        %%         ')
    print('  % %%%%%%%   %%%%%   %%%%%%%%%%%%    %%%    %%%%%%%%      ')
    print('  ---------------------------------------------------------')
    print('                Simple MOC Ray-Tracing code                ')
    print('                     Version: 0.1-dev')
    print('                   Copyright: ZouHang')
    print('')
