function [decoder_output,survivor_state,cumulated_metric]=viterbi_s(G,k,channel_output)
%VITERBI	The Viterbi decoder for convolutional codes
%		[decoder_output,survivor_state,cumulated_metric]=viterbi(G,k,channel_output)
%		G is a n x Lk matrix each row of which
%          	determines the connections from the shift register to the
%          	n-th output of the code, k/n is the rate of the code.
%          	survivor_state is a matrix showing the optimal path through
%          	the trellis. The metric is given in a separate function metric(x,y)
%          	and can be specified to accomodate hard and soft decision.
%          	This algorithm minimizes the metric rather than maximizing
%          	the likelihood.
 
    n=size(G,1);
    %  check the sizes
    if rem(size(G,2),k) ~=0 
      error('Size of G and k do not agree')
    end
    if rem(size(channel_output,2),n) ~=0
      error('channel output not of the right size')
    end
    L=size(G,2)/k;
    number_of_states=2^((L-1)*k);
    %  generate state transition matrix, output matrix, and input matrix
    for j=0:number_of_states-1
      for l=0:2^k-1
        [next_state,memory_contents]=nxt_stat(j,l,L,k);
        input(j+1,next_state+1)=l;
        branch_output=rem(memory_contents*G',2);
        nextstate(j+1,l+1)=next_state;
        output(j+1,l+1)=bin2deci(branch_output);
      end
    end
    state_metric=zeros(number_of_states,2);
    depth_of_trellis=length(channel_output)/n;
    channel_output_matrix=reshape(channel_output,n,depth_of_trellis);
    survivor_state=zeros(number_of_states,depth_of_trellis+1);
    %  start decoding of non-tail channel outputs
    for i=1:depth_of_trellis-L+1
      flag=zeros(1,number_of_states);
      if i <= L
        step=2^((L-i)*k);
      else
        step=1;
      end
      for j=0:step:number_of_states-1
        for l=0:2^k-1
          branch_metric=0;
          binary_output=deci2bin(output(j+1,l+1),n);
          for ll=1:n
            branch_metric=branch_metric+metric_s(channel_output_matrix(ll,i),2*binary_output(ll)-1);
          end
          if((state_metric(nextstate(j+1,l+1)+1,2) > state_metric(j+1,1)...
            +branch_metric) | flag(nextstate(j+1,l+1)+1)==0)
            state_metric(nextstate(j+1,l+1)+1,2) = state_metric(j+1,1)+branch_metric;
            survivor_state(nextstate(j+1,l+1)+1,i+1)=j;
            flag(nextstate(j+1,l+1)+1)=1;
          end
        end
      end
      state_metric=state_metric(:,2:-1:1);
    end
    %  start decoding of the tail channel-outputs
    for i=depth_of_trellis-L+2:depth_of_trellis
      flag=zeros(1,number_of_states);
      last_stop=number_of_states/(2^((i-depth_of_trellis+L-2)*k));
      for j=0:last_stop-1
          branch_metric=0;
          binary_output=deci2bin(output(j+1,1),n);
          for ll=1:n
            branch_metric=branch_metric+metric_s(channel_output_matrix(ll,i),2*binary_output(ll)-1);
          end
          if((state_metric(nextstate(j+1,1)+1,2) > state_metric(j+1,1)...
            +branch_metric) | flag(nextstate(j+1,1)+1)==0)
            state_metric(nextstate(j+1,1)+1,2) = state_metric(j+1,1)+branch_metric;
            survivor_state(nextstate(j+1,1)+1,i+1)=j;
            flag(nextstate(j+1,1)+1)=1;
          end
      end
      state_metric=state_metric(:,2:-1:1);
    end
    %  generate the decoder output from the optimal path
    state_sequence=zeros(1,depth_of_trellis+1);
    state_sequence(1,depth_of_trellis)=survivor_state(1,depth_of_trellis+1);
    for i=1:depth_of_trellis
      state_sequence(1,depth_of_trellis-i+1)=survivor_state((state_sequence(1,depth_of_trellis+2-i)...
      +1),depth_of_trellis-i+2);
    end
    decodeder_output_matrix=zeros(k,depth_of_trellis-L+1);
    for i=1:depth_of_trellis-L+1
      dec_output_deci=input(state_sequence(1,i)+1,state_sequence(1,i+1)+1);
      dec_output_bin=deci2bin(dec_output_deci,k);
      decoder_output_matrix(:,i)=dec_output_bin(k:-1:1)';
    end
    decoder_output=reshape(decoder_output_matrix,1,k*(depth_of_trellis-L+1));
    cumulated_metric=state_metric(1,1);
end