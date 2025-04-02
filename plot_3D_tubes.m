tic
figure(3)
hold on
grid on

for i=1:length(Vessels)
    for j=1:size(Vessels{1,i},2)
        
        D=Vessels{1,i}{1,j}(1,:)-Vessels{1,i}{1,j}(2,:);
        R=Vessels{1,i}{2,j}(1,1);
        [X, Y, Z] = cylinder2P([R R], 5, 2, Vessels{1,i}{1,j}(1,:),Vessels{1,i}{1,j}(2,:));
        C=(56/255)*ones(2,5);
        C(:,:,2)=(75/255)*ones(2,5);
        C(:,:,3)=(137/255)*ones(2,5);
        surf(X, Y, Z, C);
        grid off
    end
end

view(0,0)
%print('st.eps','-depsc')

toc     
        
        
        
        
