figure
for i = 1:36
	subplot(6,6,i)
	plot(rf,i);
end

subplot(6,6,1)
ylabel('Frequency')
subplot(6,6,7)
ylabel('Frequency')
subplot(6,6,13)
ylabel('Frequency')
subplot(6,6,19)
ylabel('Frequency')
subplot(6,6,25)
ylabel('Frequency')
subplot(6,6,31)
ylabel('Frequency')
xlabel('ISI (ms)')
subplot(6,6,32)
xlabel('ISI (ms)')
subplot(6,6,33)
xlabel('ISI (ms)')
subplot(6,6,34)
xlabel('ISI (ms)')
subplot(6,6,35)
xlabel('ISI (ms)')
subplot(6,6,36)
xlabel('ISI (ms)')

% arrange ISI's according to yscale
% cells with counts over 100
figure
cell_order = [2; 5; 6; 8; 10; 11; 12; 13; 14; 15; 17; 20; 23; 28; 29; 30; 31; 32; 34; 35; 36];
for i = 1:21
	subplot(6,7,i)
	plot(rf,cell_order(i),'xMax',15);
end

% cells with counts over 10
cell_order2 = [1; 3; 4; 7; 9; 16; 18; 21; 22; 24; 25; 26; 27; 33];
for i = 1:14
	subplot(6,7,i+21)
	plot(rf,cell_order2(i),'xMax',15);
end

% cells with counts under 10
cell_order3 = 19;
subplot(6,7,36)
plot(rf,cell_order3(1),'xMax',15);


% plot recovery functions below ISIs
% cells with counts over 100
figure
for i = 1:21
	n = floor((i-1)/7);
	subplot(6,7,n*7+i)
	plot(rf,cell_order(i),'xMax',15);
	subplot(6,7,(n+1)*7+i)
	plot(rf,cell_order(i),'xMax',15,'recovery');
end
subplot(6,7,1)
ylabel('Frequency')
subplot(6,7,8)
ylabel('Recovery')
subplot(6,7,15)
ylabel('Frequency')
subplot(6,7,22)
ylabel('Recovery')
subplot(6,7,29)
ylabel('Frequency')
subplot(6,7,36)
ylabel('Recovery')
xlabel('ISI (ms)')

% cells with counts over 10
figure
for i = 1:14
	n = floor((i-1)/7);
	subplot(6,7,n*7+i)
	plot(rf,cell_order2(i),'xMax',15);
	subplot(6,7,(n+1)*7+i)
	plot(rf,cell_order2(i),'xMax',15,'recovery');
end

% cells with counts under 10
%figure
subplot(6,7,29)
plot(rf,cell_order3(1),'xMax',15);
% [rec,edges] = relrefrac(jonscells,cell_order3(1));
subplot(6,7,36)
plot(rf,cell_order3(1),'xMax',15,'recovery');
% plot(edges(1:length(rec)),rec,'d-')
% xlim([0 3])
subplot(6,7,1)
ylabel('Frequency')
subplot(6,7,8)
ylabel('Recovery')
subplot(6,7,15)
ylabel('Frequency')
subplot(6,7,22)
ylabel('Recovery')
subplot(6,7,29)
ylabel('Frequency')
subplot(6,7,36)
ylabel('Recovery')
xlabel('ISI (ms)')

% get surrogate spikes
qt = rfd.qt(:,1);
wt = rfd.wt{1};
% set number of points to step
cstep = length(wt);
% create sptrain so we don't have to keep changing memory size
% assume you can't have more than 1 spike in 1 ms so make sptrain
% equal to the duration of a repetition in ms
sptrain = zeros(rfd.duration,1);
% use spt to check to see if we are at the end of one repetition
spt = 0;
spti = 1;
% calculate running sum for cstep points
matSum = tril(ones(cstep,cstep));

% get first spike by using w(t) = 1
% get random number
r = rand;
rln = -log(r);
% initialize cend for loop 
cend = 0;
eValue = 0;
cSi = [];
while(1)
	while( isempty(cSi) )
		% get the next start and end values
		cstart = cend + 1;
		cend = cend + cstep;
		% compute the next cstep cummulative sums
		cSum = eValue + matSum * qt(cstart:cend);
		% find if there is an index greater than rln
		cSi = find(cSum>rln);
		% get the last value from the cummulative sum for next
		% calculation
		eValue = cSum(cstep);
	end
	% value found so figure out where to put the spike
	% get index that first exceeds rln
	% if it was index 22, cstart will be 21, cSi(1) will be 2 so in 
	% order to get back 22, we subtract 1 from cstart + cSi(1)
	xrlni = cstart+cSi(1)-1;
	spt = rfd.rtEdges(xrlni);
	if( spt < rfd.duration )
		% add value to spike train
		sptrain(spti) = spt;
		% increment spti
		spti = spti + 1;
		% get new random number
		r = rand;
		rln =  -log(r);
		% reset cstart, cend, eValue and cSi
		cstart = xrlni;
		% if cstep is 10 and cstart is 22 then cend = 22 + 10 - 1
		% will be 10 values
		cend = cstart + cstep - 1;
		% find cummulative sum again taking into account the
		% relative refractory period
		cSum = matSum * (qt(cstart:cend) .* wt);
		% find if there is an index greater than rln
		cSi = find(cSum>rln);
		% get the last value from the cummulative sum for next
		% calculation
		eValue = cSum(cstep);
	else
		% break out of while loop
		break;
	end
end	% end of while(1) loop
