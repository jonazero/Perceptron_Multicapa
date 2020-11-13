classdef PMulticapa_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        UIAxes                        matlab.ui.control.UIAxes
        ClasificarButton              matlab.ui.control.Button
        TasadeaprendizajeSliderLabel  matlab.ui.control.Label
        TasadeaprendizajeSlider       matlab.ui.control.Slider
        BorrarButton                  matlab.ui.control.Button
        ClaseButtonGroup              matlab.ui.container.ButtonGroup
        AmarilloButton                matlab.ui.control.RadioButton
        AzulButton                    matlab.ui.control.RadioButton
        NumerodeneuronasenlacapaocultaSpinnerLabel  matlab.ui.control.Label
        NumerodeneuronasenlacapaocultaSpinner  matlab.ui.control.Spinner
        Label                         matlab.ui.control.Label
        CheckBox                      matlab.ui.control.CheckBox
    end

    properties (Access = private)
        X % Patrones
        T %Clases
    end
    methods (Access = private)    
        function Main(app)
            xrange = [0 1];
            yrange = [0 1];
            [xx1,xx2] = meshgrid(xrange(1):0.01:xrange(2), yrange(1):0.01:yrange(2));
            image_size = size(xx1);
            gridXY = [xx1(:) xx2(:)];
            BackPropagation(app,app.X,app.T, ...
                app.NumerodeneuronasenlacapaocultaSpinner.Value, ...
                xrange,yrange,image_size,gridXY);
        end
        function [Wij,Wjk] = BackPropagation(app, X, T, H,xrange,yrange,image_size,gridXY)
            
            %% Inicializacion
            %X = Patrones
            %H = numero de neronas en la capa oculta
            X = X';
            epochs = 200000; %Numero de epocas
            errmin = 1e-9; %Error minimo deseado
            n = app.TasadeaprendizajeSlider.Value; %Tasa de aprendizaje
            [D,N] = size(X);  %D = Dimensiones de los patrones, N = numero de patrones
            X = [ones(1,N);X]; %Agrega una fila de 1 a la matriz de patrones (Umbral)
            Wij = rand(H, D+1); %Pesos aleatorios para la capa oculta
            Wjk = rand(1,H); %Pesos aleatorios para la capa de salida
            for i=1:epochs
                %%Propagacion.
                Xj = Wij*X; % Entradas por pesos de la capa oculta
                yj = Sigmoidal(app,Xj); % Se evalua la funcion sigmoidal a cada salida de la capa oculta
                Xk = Wjk*yj; %Entradas por pesos de la capa de salida
                zk = Sigmoidal(app,Xk);            
                E = T-zk; %Calculo del error
                J = mean(E.^2); %Error cuadratico medio
                %%Retropropagacion
                dk = zk.*(1-zk).*E; %ÿk
                Wjk = Wjk + n* dk*yj'; %Se actualizan los pesos de la capa de salida
                dj = (Wjk'*dk).*yj.*(1-yj); %ÿj
                Wij = Wij + n*dj*X'; %Se actualizan los pesos de la capa oculta
                if ~(mod(i,500))
                    PlotBound(app, Wij, Wjk,xrange,yrange,image_size,gridXY);
                    app.Label.Text = string(i);                
                end     
                if J < errmin
                    return
                end
            end   
        end     
        function y = Sigmoidal(~,x)
            y = 1./(1+exp(-x));
        end
        function y = Tanh(~,x)
            y = sinh(x)./(cosh(x));
        end
        function PlotLines(app,Wij)
            xx = 0:0.01:1;
            tam = length(Wij)-1;
            for i=1:tam
                plot(app.UIAxes, xx, -(Wij(i,3)/Wij(i,2))*xx-(Wij(i,1)/Wij(i,2)));
            end
        end
        function PlotBoundary(app,Wij, Wjk)
            xrange = [0 1];
            yrange = [0 1];
            [xx1,xx2] = meshgrid(xrange(1):0.01:xrange(2), yrange(1):0.01:yrange(2));
            image_size = size(xx1);
            gridXY = [xx1(:) xx2(:)];
            salida = Propagation(app,gridXY,Wij,Wjk,2);
            %salida(salida>0.5)=(1);
            %salida(salida<0.5)=(0);
            decisionmap = reshape(salida,image_size);
            imagesc(app.UIAxes, xrange, yrange, decisionmap);
            tam = length(app.X);
            for i=1:tam
                if app.T(i)==1
                   plot(app.UIAxes,app.X(i,1),app.X(i,2),'.',"Color", '#ffce00',"MarkerSize",15); 
                else
                   plot(app.UIAxes,app.X(i,1),app.X(i,2),'.',"Color", '#1979a9',"MarkerSize",15); 
                end            
            end
        end
        function PlotBound(app,Wij, Wjk,xrange,yrange,image_size,gridXY)
            salida = Propagation(app,gridXY,Wij,Wjk,2);
            %salida(salida>0.5)=(1);
            %salida(salida<0.5)=(0);
            decisionmap = reshape(salida,image_size);                     
            imagesc(app.UIAxes, xrange, yrange, decisionmap);        
            tam = length(app.X);
            for i=1:tam
                if app.T(i)==1
                   plot(app.UIAxes,app.X(i,1),app.X(i,2),'.',"Color", '#ffce00',"MarkerSize",15); 
                else
                   plot(app.UIAxes,app.X(i,1),app.X(i,2),'.',"Color", '#1979a9',"MarkerSize",15); 
                end            
            end
            pause(0.0001);  
            app.UIAxes.cla;                
        end
        function zk = Propagation(app,X,Wij,Wjk,~)
            %% Inicializacion
            %X = Patrones
            %H = numero de neronas en la capa oculta
            X = X';
            [~,N] = size(X);  %D = Dimensiones de los patrones, N = numero de patrones
            X = [ones(1,N);X]; %Agrega una fila de 1 a la matriz de patrones (Umbral)        
            %%Propagacion.
            Xj = Wij*X; % Entradas por pesos de la capa oculta
            yj = Sigmoidal(app,Xj); % Se evalua la funcion sigmoidal a cada salida de la capa oculta
            Xk = Wjk*yj; %Entradas por pesos de la capa de salida
            zk = Sigmoidal(app,Xk); %Evaluacion de la funcion Sigmoidal para la capa de salida    
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Window button down function: UIFigure
        function UIFigureWindowButtonDown(app, event)
            axis(app.UIAxes,[0 1 0 1]);
            temp = app.UIAxes.CurrentPoint; % Returns 2x3 array of points
            loc = [temp(1,1) temp(1,2)]; % Gets the (x,y) coordinates
            hold(app.UIAxes, 'on' )
            if app.AmarilloButton.Value == true
                plot(app.UIAxes,loc(1),loc(2),'.',"Color", '#f8f916',"MarkerSize",15);
                app.X = [app.X; loc(1),loc(2)];
                app.T = [app.T, 1];
            elseif app.AzulButton.Value == true
                plot(app.UIAxes,loc(1),loc(2),'.',"Color", '#3d26a8',"MarkerSize",15);
                app.X = [app.X; loc(1),loc(2)];
                app.T = [app.T, 0];
            end
            app.ClasificarButton.Enable = true;
        end

        % Button pushed function: ClasificarButton
        function ClasificarButtonPushed(app, event)
            Main(app)
        end

        % Button pushed function: BorrarButton
        function BorrarButtonPushed(app, event)
            app.X = [];
            app.T = [];
            app.UIAxes.cla;
            app.ClasificarButton.Enable = false;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 569 706];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.WindowButtonDownFcn = createCallbackFcn(app, @UIFigureWindowButtonDown, true);
            app.UIFigure.Scrollable = 'on';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Red neuronal multicapa')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.PlotBoxAspectRatio = [1.1508120649652 1 1];
            app.UIAxes.Position = [1 194 543 498];

            % Create ClasificarButton
            app.ClasificarButton = uibutton(app.UIFigure, 'push');
            app.ClasificarButton.ButtonPushedFcn = createCallbackFcn(app, @ClasificarButtonPushed, true);
            app.ClasificarButton.Enable = 'off';
            app.ClasificarButton.Position = [41 78 121 56];
            app.ClasificarButton.Text = 'Clasificar';

            % Create TasadeaprendizajeSliderLabel
            app.TasadeaprendizajeSliderLabel = uilabel(app.UIFigure);
            app.TasadeaprendizajeSliderLabel.HorizontalAlignment = 'right';
            app.TasadeaprendizajeSliderLabel.Position = [41 173 113 22];
            app.TasadeaprendizajeSliderLabel.Text = 'Tasa de aprendizaje';

            % Create TasadeaprendizajeSlider
            app.TasadeaprendizajeSlider = uislider(app.UIFigure);
            app.TasadeaprendizajeSlider.Limits = [0 0.1];
            app.TasadeaprendizajeSlider.Position = [175 182 360 3];
            app.TasadeaprendizajeSlider.Value = 0.03;

            % Create BorrarButton
            app.BorrarButton = uibutton(app.UIFigure, 'push');
            app.BorrarButton.ButtonPushedFcn = createCallbackFcn(app, @BorrarButtonPushed, true);
            app.BorrarButton.Position = [335 78 209 22];
            app.BorrarButton.Text = 'Borrar';

            % Create ClaseButtonGroup
            app.ClaseButtonGroup = uibuttongroup(app.UIFigure);
            app.ClaseButtonGroup.TitlePosition = 'centertop';
            app.ClaseButtonGroup.Title = 'Clase';
            app.ClaseButtonGroup.Position = [175 78 147 56];

            % Create AmarilloButton
            app.AmarilloButton = uiradiobutton(app.ClaseButtonGroup);
            app.AmarilloButton.Text = 'Amarillo';
            app.AmarilloButton.Position = [11 10 65 22];
            app.AmarilloButton.Value = true;

            % Create AzulButton
            app.AzulButton = uiradiobutton(app.ClaseButtonGroup);
            app.AzulButton.Text = 'Azul';
            app.AzulButton.Position = [85 10 52 22];

            % Create NumerodeneuronasenlacapaocultaSpinnerLabel
            app.NumerodeneuronasenlacapaocultaSpinnerLabel = uilabel(app.UIFigure);
            app.NumerodeneuronasenlacapaocultaSpinnerLabel.Position = [335 106 125 28];
            app.NumerodeneuronasenlacapaocultaSpinnerLabel.Text = {'Numero de neuronas'; 'en la capa oculta'};

            % Create NumerodeneuronasenlacapaocultaSpinner
            app.NumerodeneuronasenlacapaocultaSpinner = uispinner(app.UIFigure);
            app.NumerodeneuronasenlacapaocultaSpinner.Limits = [2 10];
            app.NumerodeneuronasenlacapaocultaSpinner.HorizontalAlignment = 'center';
            app.NumerodeneuronasenlacapaocultaSpinner.Position = [450 106 94 28];
            app.NumerodeneuronasenlacapaocultaSpinner.Value = 4;

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.Position = [42 29 82 22];

            % Create CheckBox
            app.CheckBox = uicheckbox(app.UIFigure);
            app.CheckBox.Position = [251 29 80 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = PMulticapa_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end